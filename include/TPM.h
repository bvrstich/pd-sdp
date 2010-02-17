#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;

class TPM : public Matrix {

   friend ostream &operator<<(ostream &,TPM &);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(TPM &);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(TPM &b,SUP &D);

      //los het stelsel op
      int solve(TPM &b,SUP &D);

      //de G down en inverse G up
      void G(int option,PHM &);

      void min_unit(double scale);

      void min_qunit(double scale);

   private:

      static int **t2s;
      static int **s2t;

      static int counter;

      int N;//nr of particles
      int M;//dim sp space
      int n;//dim tp space

};

#endif
