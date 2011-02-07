#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M
 * @param M dimension of single particle space and dimension of the Matrix
 * @param N Nr of particles
 */
GutMat::GutMat(int M,int N) : Matrix(M) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param gm_copy content of this matrix will be copied into the constructed matrix
 */
GutMat::GutMat(const GutMat &gm_copy) : Matrix(gm_copy) {

   this->M = gm_copy.gM();
   this->N = gm_copy.gN();

}

/**
 * destructor
 */
GutMat::~GutMat(){

}

/**
 * @return nr of particles
 */
int GutMat::gN() const
{
   return N;
}

/**
 * @return dimension of sp space
 */
int GutMat::gM() const
{
   return M;
}

ostream &operator<<(ostream &output,GutMat &gm_p){

   for(int i = 0;i < gm_p.M;++i)
      for(int j = 0;j < gm_p.M;++j)
         output << i << "\t" << j << "\t" << gm_p(i,j) << endl;

   return output;

}

/**
 * makes the approximated spm on singly occupied space
 * @param tpm The input TPM
 */
void GutMat::p(const TPM &tpm){

   *this = 0.0;

   SPM spm(1.0/(N - 1.0),tpm);

   for(int a = 0;a < M;++a){

      (*this)(a,a) = spm(Math::adjoint(a),Math::adjoint(a));

      for(int b = a;b < M;++b)
         (*this)(a,b) += spm(a,b) - tpm(a,Math::adjoint(b),b,Math::adjoint(b)) - tpm(a,Math::adjoint(a),b,Math::adjoint(a));

   }

   this->symmetrize();

}

/**
 * makes the approximated one hole matrix on singly occupied space
 * @param tpm The input TPM
 */
void GutMat::q(const TPM &tpm){

   *this = 0.0;

   SPM spm(1.0/(N - 1.0),tpm);

   double ward = 2.0*tpm.trace()/(N*(N - 1.0));

   for(int a = 0;a < M;++a){

      (*this)(a,a) = ward - spm(Math::adjoint(a),Math::adjoint(a));

      for(int b = a;b < M;++b)
         (*this)(a,b) += tpm(Math::adjoint(b),b,Math::adjoint(b),a) + tpm(Math::adjoint(a),b,Math::adjoint(a),a) - spm(a,b);

   }

   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
