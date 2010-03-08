#include <iostream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;

#include "include.h"

/**
 * constructor 
 * @param n dimension of the matrix
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
 * @param mat_copy The matrix you want to be copied into the object you are constructing
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
 * overload the equality operator
 * @param matrix_copy The matrix you want to be copied into this
 */
Matrix &Matrix::operator=(Matrix &matrix_copy){

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,matrix_copy.matrix[0],&incx,matrix[0],&incy);

   return *this;

}

/**
 * Make all the number in your matrix equal to the number a, e.g. usefull for initialization (Matrix M = 0)
 * @param a the number
 */
Matrix &Matrix::operator=(double a){

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         matrix[j][i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param matrix_pl The matrix you want to add to this
 */
Matrix &Matrix::operator+=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param matrix_pl The matrix you want to deduct from this
 */
Matrix &Matrix::operator-=(Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param matrix_pl the Matrix to be multiplied by alpha and added to this
 */
Matrix &Matrix::daxpy(double alpha,Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}
/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
Matrix &Matrix::operator/=(double c){

   int dim = n*n;
   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&dim,&alpha,matrix[0],&inc);

   return *this;

}

/**
 * write access to your matrix, change the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double &Matrix::operator()(int i,int j){

   return matrix[j][i];

}

/**
 * read access to your matrix, view the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double Matrix::operator()(int i,int j) const {

   return matrix[j][i];

}

/**
 * @return the underlying pointer to matrix, useful for mkl applications
 */
double **Matrix::gMatrix(){

   return matrix;

}

/**
 * @return the dimension of the matrix
 */
int Matrix::gn(){

   return n;

}

/**
 * @return the trace of the matrix:
 */
double Matrix::trace(){

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += matrix[i][i];

   return ward;

}

/**
 * Diagonalizes symmetric matrices. Watch out! The current matrix (*this) is destroyed, in it
 * the eigenvectors will be stored (one in every column).
 * @param eigenvalues the pointer of doubles in which the eigenvalues will be storen, Watch out, its memory
 * has to be allocated on the dimension of the matrix before you call the function.
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
 * @return inproduct of (*this) matrix with matrix_i, defined as Tr (A B)
 * @param matrix_i input matrix
 */
double Matrix::ddot(Matrix &matrix_i){

   int dim = n*n;
   int inc = 1;

   return ddot_(&dim,matrix[0],&inc,matrix_i.matrix[0],&inc);

}

/**
 * Invert positive semidefinite symmetric matrix is stored in (*this), original matrix (*this) is destroyed
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
 * Scale the matrix (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Matrix::dscal(double alpha){

   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&alpha,matrix[0],&inc);

}

/**
 * Fill the matrix with random numbers.
 */
void Matrix::fill_Random(){

   srand(time(NULL));

   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         matrix[j][i] = (double) rand()/RAND_MAX;

   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[i][j] = matrix[j][i];

}

/**
 * Take the square root out of the positive semidefinite matrix, destroys original matrix, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
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
 * Multiply this matrix with diagonal matrix
 * @param diag Diagonal matrix to multiply with this, has to be allocated on matrix dimension.
 */
void Matrix::mdiag(double *diag){

   int inc = 1;

   for(int i = 0;i < n;++i)
      dscal_(&n,diag + i,matrix[i],&inc);

}

/**
 * Multiply symmetric matrix object left en right with symmetric matrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map matrix that will be multiplied to the left en to the right of matrix object
 * @param object central matrix
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
 * Matrix product of two general matrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
Matrix &Matrix::mprod(Matrix &A, Matrix &B){

   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&trans,&trans,&n,&n,&n,&alpha,A.matrix[0],&n,B.matrix[0],&n,&beta,matrix[0],&n);

   return *this;

}

/**
 * Copy upper triangle into lower triangle.
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