#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "SPM.h"
#include "lapack.h"

/**
 * constructor, maakt een Matrix object aan met dimensie M
 * @param M aantal sp orbitals, tevens dimensie van de matrix
 * @param N aantal deeltjes
 */
SPM::SPM(int M,int N) : Matrix(M) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor, maakt een Matrix object aan dat een kopie is van de inputmatrix
 * @param spm_copy inputmatrix
 */
SPM::SPM(SPM &spm_copy) : Matrix(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

//!destructor
SPM::~SPM(){

}

/**
 * @return aantal deeltjes
 */

int SPM::gN(){

   return N;

}

/**
 * @return aantal sp orbitals en dimensie van de matrix
 */

int SPM::gM(){

   return M;

}

ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.M;++i)
      for(int j = 0;j < spm_p.M;++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}
