#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * Dit is een zelfgeschreven matrix-klasse voor vierkante matrices. Het is een wrapper rond een dubbele pointer
 * met veelgebruike lapack en blas routines als voorgeprogrammeerde memberfuncties.
 */
class Matrix{

   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << matrix << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << matrix << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param matrix_p de te printen matrix
    */
   friend ostream &operator<<(ostream &output,Matrix &matrix_p);

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

      //!De dubbele pointer van double's die de eigenlijke matrix bevat.
      double **matrix;

      //!De dimensie van de vierkante matrix
      int n;

};

#endif
