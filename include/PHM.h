#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"
#include "SPM.h"

class PHM : public Matrix {

   friend ostream &operator<<(ostream &,PHM &);

   public:
      
      //constructor
      PHM(int M,int N);

      //copy constructor
      PHM(PHM &);

      //destructor
      virtual ~PHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      //double operator()(int a,int b,int c,int d) const;

      //change the numbers in sp mode
      double &operator()(int a,int b,int c,int d);

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      //geef dim terug
      int gn();

      void G(int option,TPM &);

      double skew_trace();

      void min_gunit(double scale);

   private:

      //counter of objects
      static int counter;

      //relationship bewtween ph and sp basis
      static int **ph2s;
      static int **s2ph;

      int N;//nr of particles
      int M;//dim sp space
      int n;//dim ph space

};

#endif
