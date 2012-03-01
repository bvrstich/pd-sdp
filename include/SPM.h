#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This class SPM was written for single particle matrices. It inherits from the class Matrix and expands it with
 * specific memberfunction and a knowledge of the nr of sp orbitals and particles.
 */

class SPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //destructor
      virtual ~SPM();

      using Matrix::operator=;

      int gN() const;

      int gM() const;

      /**
       * constructs a SPM from a TPM or a PHM. Definition of these functions are in different notes.
       * Actually this function is defined as bar * scale.
       * @param scale The factor that multiplies the bar(MT), e.g. 1/(N - 1) for a normal single particle density matrix
       * @param MT PHM or TPM inputmatrix.
       */
      template<class MatrixType>
         void bar(double scale,const MatrixType &MT){

            for(int a = 0;a < MT.gM();++a)
               for(int b = a;b < MT.gM();++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < MT.gM();++l)
                     (*this)(a,b) += MT(a,l,b,l);

                  (*this)(a,b) *= scale;

               }

            this->symmetrize();

         }
      void bar(const T2PM &MT);

      static void init(int,int);

   private:

      //!dimension of single particle space
      static int M;

      //!nr of particles
      static int N;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
