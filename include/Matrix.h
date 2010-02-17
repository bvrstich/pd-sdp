#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class Matrix{

   friend ostream &operator<<(ostream &,Matrix &);

   public:

      //constructor
      Matrix(int n);

      //copy constructor
      Matrix(Matrix &);

      //destructor
      virtual ~Matrix();

      //overload equality operator
      Matrix &operator=(Matrix &);

      Matrix &operator=(double );

      //overload += operator
      Matrix &operator+=(Matrix &);

      //overload -= operator
      Matrix &operator-=(Matrix &);

      Matrix &daxpy(double alpha,Matrix &);

      Matrix &mprod(Matrix &A,Matrix &B);

      Matrix &operator/=(double );

      //easy to change the numbers
      double &operator()(int i,int j);

      //easy to access the numbers
      double operator()(int i,int j) const;

      //get the pointer to the matrix
      double **gMatrix();

      int gn();

      double trace();

      //diagonaliseer symmetrische matrices
      void diagonalize(double *eigenvalues);

      double ddot(Matrix &);

      void invert();

      void dscal(double alpha);

      void fill_Random();

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void mdiag(double *diag);

      void L_map(Matrix &,Matrix &);

      void symmetrize();

   private:

      //velden
      double **matrix;

      int n;//dim of matrix

};

#endif
