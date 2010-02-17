#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "SPM.h"
#include "lapack.h"

//constructor:
SPM::SPM(int M,int N) : Matrix(M) {

   this->M = M;
   this->N = N;

}

//copy constructor()
SPM::SPM(SPM &spm_copy) : Matrix(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

//destructor
SPM::~SPM(){

}

int SPM::gN(){

   return N;

}

int SPM::gM(){

   return M;

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.M;++i)
      for(int j = 0;j < spm_p.M;++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}
