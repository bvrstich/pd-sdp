#ifndef PHPM_H
#define PHPM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 09-05-2011\n\n
 * This class, PHPM, is a class written for particle-hole-particle matrices. It is written specially for the T_3 condition. 
 * It inherits all the functions from its mother class Matrix, some special member functions and two lists that give the relationship
 * between the phh and the sp basis.
 */
class PHPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param phpm_p the PHPM you want to print
    */
   friend ostream &operator<<(ostream &output,PHPM &phpm_p);

   public:
      
      //constructor
      PHPM(int M,int N);

      //copy constructor
      PHPM(const PHPM &);

      //destructor
      virtual ~PHPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef dim terug
      int gn() const;

      //maak een PHPM van een TPM via de T3 conditie
      void T(const TPM &);

   private:

      //!static counter that counts the number of PHPM objects running in the program
      static int counter;

      //!static list of dimension [n_php][3] that takes in a php index i and returns three sp indices: a = php2s[i][0], b = php2s[i][1] and c = php2s[i][2]
      static int **php2s;

      //!static list of dimension [M][M][M] that takes three sp indices a,b and c and returns a php index i: i = s2php[a][b][c]
      static int ***s2php;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

      //!dimension of php space
      int n;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
