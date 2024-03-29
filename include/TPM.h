#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;
class DPM;
class PPHM;
class T2PM;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class TPM is a class written for two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class TPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,TPM &tpm_p);

   public:
      
      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef N terug
      int gM() const;

      void hubbard_1D(int option,double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,const TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(const TPM &b,const SUP &D);

      //los het stelsel op
      int solve(TPM &b,const SUP &D);

      //de G down en inverse G up
      void G(int option,const PHM &);

      //trace one pair of indices of DPM
      void bar(const DPM &);

      //trace one pair of indices of PPHM
      void bar(const PPHM &);

      //T1 down
      void T(int option,const DPM &);

      //T2 down
      void T(const PPHM &);

      void collaps(int option,const SUP &);

      void out(const char *filename);

      void sp_pairing(double );

      void pairing(double x[]);

      void in_sp(const char *);

      void bar(const T2PM &);

      void T(const T2PM &);

      static void init_overlap(int,int);

      static void init(int,int);

      static void clear();

   private:

      //!static list of dimension [n_tp][2] that takes in a tp index i and returns two sp indices: a = t2s[i][0] and b = t2s[i][1]
      static vector< vector<int> > t2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and returns a tp index i: i = s2t[a][b]
      static int **s2t;

      //!static variables of the inverse overlapmatrix.
      static double Sa,Sb,Sc;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
