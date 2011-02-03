#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M/2
 * @param M dimension of single particle space and dimension of the Matrix
 * @param N Nr of particles
 */
GutMat::GutMat(int M,int N) : Matrix(M/2) {

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

/* vim: set ts=3 sw=3 expandtab :*/
