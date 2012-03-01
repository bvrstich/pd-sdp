#ifndef PPHM_H
#define PPHM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 18-03-2010\n\n
 * This class, PPHM, is a class written for two-particle-one-hole matrices. It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class Matrix, some special member functions and two lists that give the relationship between the pph (two-particle one hole) and the sp basis.
 */
class PPHM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << pphm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << pphm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,PPHM &pphm_p);

   public:
      
      //constructor
      PPHM();

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //maak een PPHM van een TPM via de T2 conditie
      void T(int option,const TPM &);

      static void init(int,int);

      static void clear();

   private:

      //!static list of dimension [n_pph][3] that takes in a pph index i and returns three sp indices: a = pph2s[i][0], b = pph2s[i][1] and c = pph2s[i][2]
      static vector< vector<int> > pph2s;

      //!static list of dimension [M][M][M] that takes three sp indices a,b and c and returns a pph index i: i = s2pph[a][b][c]
      static int ***s2pph;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
