#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

int SPM::M;
int SPM::N;

/**
 * function that initializes the statics
 * @param M_i dimension of single particle space
 * @param N_i number of particles
 */
void SPM::init(int M_i,int N_i){

   M = M_i;
   N = N_i;

}

/**
 * constructor, makes matrix of dimension M
 */
SPM::SPM() : Matrix(M) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Matrix(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN() const
{
   return N;
}

/**
 * @return dimension of sp space and of matrix
 */
int SPM::gM() const
{
   return M;
}

ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.gM();++i)
      for(int j = 0;j < spm_p.gM();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

void SPM::bar(const T2PM &MT){

   for(int a = 0;a < M;a++)
      for(int b = a;b < M;b++){

         (*this)(a,b) = 0.0;

         for(int l=0;l<M;l++)
            for(int k=0;k<M;k++)
               (*this)(a,b) += MT(l,k,a,l,k,b);

         (*this)(b,a) = (*this)(a,b);

      }

}

/* vim: set ts=3 sw=3 expandtab :*/
