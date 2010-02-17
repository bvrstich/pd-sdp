#include <iostream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;

#include "lapack.h"
#include "Matrix.h"

//constructor:
Matrix::Matrix(int n){

   this->n = n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

}

//copy constructor:
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

//destructor
Matrix::~Matrix(){

   delete [] matrix[0];
   delete [] matrix;

}

//overload equality operator
Matrix &Matrix::operator=(Matrix &matrix_copy){

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,matrix_copy.matrix[0],&incx,matrix[0],&incy);

   return *this;

}

Matrix &Matrix::operator=(double a){

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         matrix[j][i] = a;

   return *this;

}

//overload += operator
Matrix &Matrix::operator+=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

//overload -= operator
Matrix &Matrix::operator-=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

//+= times a constant 
Matrix &Matrix::daxpy(double alpha,Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

// divide by a constant
Matrix &Matrix::operator/=(double c){

   int dim = n*n;
   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&dim,&alpha,matrix[0],&inc);

   return *this;

}

//change the numbers
double &Matrix::operator()(int i,int j){

   return matrix[j][i];

}

//access the numbers
double Matrix::operator()(int i,int j) const {

   return matrix[j][i];

}

//get pointer to matrix matrix
double **Matrix::gMatrix(){

   return matrix;

}

int Matrix::gn(){

   return n;

}

double Matrix::trace(){

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += matrix[i][i];

   return ward;

}

//diagonaliseer symmetrische matrices
void Matrix::diagonalize(double *eigenvalues){

   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,matrix[0],&n,eigenvalues,work,&lwork,&info);

   delete [] work;

}

double Matrix::ddot(Matrix &matrix_i){

   int dim = n*n;
   int inc = 1;

   return ddot_(&dim,matrix[0],&inc,matrix_i.matrix[0],&inc);

}

void Matrix::invert(){

   char uplo = 'U';

   int INFO;

   dpotrf_(&uplo,&n,matrix[0],&n,&INFO);//cholesky decompositie

   dpotri_(&uplo,&n,matrix[0],&n,&INFO);//inverse berekenen

   //terug symmetrisch maken:
   this->symmetrize();

}

void Matrix::dscal(double alpha){

   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&alpha,matrix[0],&inc);

}

void Matrix::fill_Random(){

   srand(time(NULL));

   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         matrix[j][i] = (double) rand()/2147483647.0;

   this->symmetrize();

}

//^{1/2} of ^{-1/2}
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

void Matrix::mdiag(double *diag){

   int inc = 1;

   for(int i = 0;i < n;++i)
      dscal_(&n,diag + i,matrix[i],&inc);

}

//symmetrische matrix links en rechts vermenigvuldigen met symmetrische matrix:
//this = map*object*map
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

Matrix &Matrix::mprod(Matrix &A,Matrix &B){

   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&trans,&trans,&n,&n,&n,&alpha,A.matrix[0],&n,B.matrix[0],&n,&beta,matrix[0],&n);

   return *this;

}

//kopieer bovendriehoek in benedendriehoek
void Matrix::symmetrize(){

   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[i][j] = matrix[j][i];

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,Matrix &matrix_p){

   for(int i = 0;i < matrix_p.gn();++i)
      for(int j = 0;j < matrix_p.gn();++j)
         output << i << "\t" << j << "\t" << matrix_p(i,j) << endl;

   return output;

}
