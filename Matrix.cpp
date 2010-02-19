#include <iostream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;

#include "lapack.h"
#include "Matrix.h"

/**
 * constructor
 * @param n dimensie van de matrix
 */

Matrix::Matrix(int n){

   this->n = n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

}

/**
 * copy constructor
 * @param mat_copy Deze matrix zal gekopieerd worden in het nieuw aangemaakte object.
 */

Matrix::Matrix(Matrix &mat_copy){

   this->n = mat_copy.n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,mat_copy.matrix[0],&incx,matrix[0],&incy);

}

/**
 * Destructor
 */

Matrix::~Matrix(){

   delete [] matrix[0];
   delete [] matrix;

}

/**
 * overload equality operator
 * @param matrix_copy Deze matrix zal gekopieerd worden in this.
 */

Matrix &Matrix::operator=(Matrix &matrix_copy){

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,matrix_copy.matrix[0],&incx,matrix[0],&incy);

   return *this;

}

/**
 * Zet alle getallen in de matrix gelijk aan a
 * @param a Het getal
 */

Matrix &Matrix::operator=(double a){

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         matrix[j][i] = a;

   return *this;

}

/**
 * Overload += operator
 * @param matrix_pl De matrix die opgeteld moet worden bij this
 */
Matrix &Matrix::operator+=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * Overload -= operator
 * @param matrix_pl De matrix die afgetrokken moet worden van this
 */

Matrix &Matrix::operator-=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * Bereken deze matrix plus een constante maal een andere matrix
 * @param alpha Het getal waarmee matrix_pl vermenigvuldig wordt
 * @param matrix_pl De matrix die alpha maal opgeteld wordt bij this
 */

Matrix &Matrix::daxpy(double alpha,Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * Overload /= operator
 * @param c Deel alle getallen in de matrix door c
 */

Matrix &Matrix::operator/=(double c){

   int dim = n*n;
   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&dim,&alpha,matrix[0],&inc);

   return *this;

}

/**
 * Overloaded () operator, write toegang tot de getallen in de matrix, om het makkelijk te maken
 * om te werken met lapack en blas routines worden deze matrices automatisch getransponeerd. We kunnen 
 * dus werken met de C notatie voor matrices terwijl lapack werkt met de fortran notatie.
 * @param i rij
 * @param j kolom 
 * @return wijzigbaar getal op de plaats (i,j) in de matrix
 */
double &Matrix::operator()(int i,int j){

   return matrix[j][i];

}

/**
 * Overloaded () operator, read toegang tot de getallen in de matrix, om het makkelijk te maken
 * om te werken met lapack en blas routines worden deze matrices automatisch getransponeerd. We kunnen 
 * dus werken met de C notatie voor matrices terwijl lapack werkt met de fortran notatie.
 * @param i rij
 * @param j kolom 
 * @return niet-wijzigbaar getal op de plaats (i,j) in de matrix
 */

double Matrix::operator()(int i,int j) const {

   return matrix[j][i];

}

/**
 * Geeft de eigenlijke dubbele pointer terug, soms handig wanneer met mkl gewerkt wordt.
 * @return dubbele pointer van double, matrix.
 */

double **Matrix::gMatrix(){

   return matrix;

}

/**
 * @return dimensie van de matrix
 */

int Matrix::gn(){

   return n;

}

/**
 * berekend de trace van de matrix
 * @return trace van de matrix
 */

double Matrix::trace(){

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += matrix[i][i];

   return ward;

}

/**
 * Diagonaliseerd symmetrische matrices. Opgelet, huidige matrix wordt vernietigd, 
 * daar zullen de eigenvectoren opgeslaan worden (In de kolommen!)
 * @param eigenvalues de pointer van doubles waarin de eigenwaarden opgeslaan gaan worden, 
 * geheugen moet op voorhand gealloceerd zijn op de matrixdimensie!
 */

void Matrix::diagonalize(double *eigenvalues){

   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,matrix[0],&n,eigenvalues,work,&lwork,&info);

   delete [] work;

}

/**
 * Berekend het inproduct van twee matrices gedefinieerd als Tr(M_1 M_2)
 * @param matrix_i De matrix waarmee het inproduct van this genomen wordt
 * @return double met inproduct in.
 */

double Matrix::ddot(Matrix &matrix_i){

   int dim = n*n;
   int inc = 1;

   return ddot_(&dim,matrix[0],&inc,matrix_i.matrix[0],&inc);

}

/**
 * Inverteer de symmetrische, positief definiete matrix.
 * Gebruikt de lapack implementatie van cholesky decompostie dus de matrix MOET
 * positief definiet zijn! Vernietigd de oorspronkelijke matrix.
 */

void Matrix::invert(){

   char uplo = 'U';

   int INFO;

   dpotrf_(&uplo,&n,matrix[0],&n,&INFO);//cholesky decompositie

   dpotri_(&uplo,&n,matrix[0],&n,&INFO);//inverse berekenen

   //terug symmetrisch maken:
   this->symmetrize();

}

/**
 * Herschaal de matrix met een factor alpha
 * @param alpha Het bewuste getal
 */
void Matrix::dscal(double alpha){

   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&alpha,matrix[0],&inc);

}

/**
 * Vul de matrix met random getallen, gebruikt de standaard randomnummergenerator van C++
 * Wordt geinitialiseerd (geseed) op de tijd dus kan twee maal zelfde resultaat reeks geven
 * als de computertijd hetzelfde is.
 */

void Matrix::fill_Random(){

   srand(time(NULL));

   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         matrix[j][i] = (double) rand()/2147483647.0;

   this->symmetrize();

}

/**
 * Neem de positieve of negatieve vierkantswortel uit de matrix. Opgelet, matrix wordt vernietigd
 * @param option = +1 Neem dan de Positieve vierkantswortel, = -1 Neem dan de negatieve vierkantswortel
 */

void Matrix::sqrt(int option){

   Matrix hulp(*this);

   double *eigen = new double [n];

   hulp.diagonalize(eigen);

   if(option == 1)
      for(int i = 0;i < n;++i)
         eigen[i] = std::sqrt(eigen[i]);
   else
      for(int i = 0;i < n;++i)
         eigen[i] = 1.0/std::sqrt(eigen[i]);

   //hulp opslaan
   Matrix hulp_c = hulp;

   //vermenigvuldigen met diagonaalmatrix
   hulp_c.mdiag(eigen);

   //en tenslotte de laatste matrixvermenigvuldiging
   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&transA,&transB,&n,&n,&n,&alpha,hulp_c.matrix[0],&n,hulp.matrix[0],&n,&beta,matrix[0],&n);

   delete [] eigen;

}

/**
 * Vermenigvuldig this met een diagonaalmatrix
 * @param diag de diagonaalmatrix waarmee this vermenigvuldig wordt
 */

void Matrix::mdiag(double *diag){

   int inc = 1;

   for(int i = 0;i < n;++i)
      dscal_(&n,diag + i,matrix[i],&inc);

}

/**
 * symmetrische matrix links en rechts vermenigvuldigen met symmetrische matrix:\n
 * this = map*object*map
 * @param map de matrix die links en rechts inwerkt op object
 * @param object wordt links en rechts vermenigvuldigd met map
 */
void Matrix::L_map(Matrix &map,Matrix &object){
   
   char side = 'L';
   char uplo = 'U';

   double alpha = 1.0;
   double beta = 0.0;

   double *hulp = new double [n*n];

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix[0],&n,object.matrix[0],&n,&beta,hulp,&n);

   side = 'R';

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix[0],&n,hulp,&n,&beta,matrix[0],&n);

   delete [] hulp;

   //expliciet symmetriseren van de uit matrix
   this->symmetrize();

}

/**
 * algemeen matrixproduct tussen twee matrices
 * 
 * @param A linkse matrix
 * @param B rechtse matrix
 * @return Het product wordt teruggegeven
 */

Matrix &Matrix::mprod(Matrix &A,Matrix &B){

   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&trans,&trans,&n,&n,&n,&alpha,A.matrix[0],&n,B.matrix[0],&n,&beta,matrix[0],&n);

   return *this;

}

/**
 * Kopieer bovendriehoek in benedendriehoek van matrix.
 */
void Matrix::symmetrize(){

   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[i][j] = matrix[j][i];

}

ostream &operator<<(ostream &output,Matrix &matrix_p){

   for(int i = 0;i < matrix_p.gn();++i)
      for(int j = 0;j < matrix_p.gn();++j)
         output << i << "\t" << j << "\t" << matrix_p(i,j) << endl;

   return output;

}
