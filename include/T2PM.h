#ifndef T2PM_H
#define T2PM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"

/**
 * @author Ward Poelmans
 * @date 11-05-2010\n\n
 */
class T2PM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << T2PM_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << T2PM_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param T2PM_p the T2PM you want to print
    */
   friend ostream &operator<<(ostream &output,const T2PM &T2PM_p);

   public:

      //constructor
      T2PM();

      //copy constructor
      T2PM(const T2PM &);

      //destructor
      virtual ~T2PM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d,int e,int f) const;

      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //maak een T2PM van een TPM
      void T(const TPM &);

      operator PPHM() const;

      operator SPM() const;

      static void init(int,int);

      static void clear();

   private:

      //!static list of dimension [n_dp][3] that takes in a pph index i and returns three sp indices: a = dp2s[i][0], b = dp2s[i][1] and c = dp2s[i][2]
      static vector< vector<int> > pph2s;

      //!static list of dimension [M][M][M] that takes three sp indices a,b and c and returns a pph index i: i = s2dp[a][b][c]
      static int ***s2pph;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif /* T2PM_H */

/* vim: set ts=3 sw=3 expandtab :*/
