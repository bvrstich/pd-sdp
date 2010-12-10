#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int LinIneq::counter = 0;

LinCon **LinIneq::li;

/**
 * constructor of a LinIneq object
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param nr nr of constraints
 */
LinIneq::LinIneq(int M,int N,int nr){

   this->M = M;
   this->N = N;

   this->nr = nr;

   if(counter == 0){

      li = new LinCon * [nr];

      for(int i = 0;i < nr;++i)
         li[i] = new LinCon(M,N);

   }

   proj = new double [nr];

   ++counter;

}

/**
 * copy constructor
 * @param li_copy The LinIneq object to be copied
 */
LinIneq::LinIneq(const LinIneq &li_copy){

   this->M = li_copy.gM();
   this->N = li_copy.gN();

   this->nr = li_copy.gnr();

   this->tr = li_copy.gtr();

   if(counter == 0){

      li = new LinCon * [nr];

      for(int i = 0;i < nr;++i)
         li[i] = new LinCon(li_copy[i]);

   }

   proj = new double [nr];

   for(int i = 0;i < nr;++i)
      proj[i] = li_copy.gproj(i);

   ++counter;

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   delete [] proj;

   if(counter == 1){

      for(int i = 0;i < nr;++i)
         delete li[i];

      delete [] li;

   }

   counter--;

}

/**
 * @return the number of currently applied linear constraints
 */
int LinIneq::gnr() const{

   return nr;

}

/**
 * read and write access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
LinCon &LinIneq::operator[](int i){

   return *li[i];

}

/**
 * read only access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
const LinCon &LinIneq::operator[](int i) const{

   return *li[i];

}

/**
 * @return nr of particles
 */
int LinIneq::gN() const{

   return N;

}

/**
 * @return nr of sp orbitals
 */
int LinIneq::gM() const{

   return M;

}

/**
 * @param li_epsi the LinIneq object filled with the Newton direction found by the cg loop.
 * @return the first singularity in the potential along the "epsilon" direction.
 */
double LinIneq::min_end(const LinIneq &li_epsi) const{

   double min_end = -constraint(0)/li_epsi.constraint(0);

   if(min_end < 0.0)
      min_end = 1.0e+15;

   for(int i = 1;i < nr;++i){

      double tmp = -constraint(i)/li_epsi.constraint(i);

      if(tmp < 0.0)
         tmp = 1.0e+15;

      if(tmp < min_end)
         min_end = tmp;

   }

   return min_end;

}

/**
 * Calculate the projection of the input TPM object on the constraint matrices of the elements of LinIneq
 * @param tpm The input TPM.
 */
void LinIneq::fill(const TPM &tpm){

   for(int i = 0;i < nr;++i)
      proj[i] = (li[i]->gI()).ddot(tpm);

   tr = 2.0*tpm.trace()/(N*(N - 1.0));

}

/**
 * @param i the index of the constraint we are interested in.
 * @return the projection on the constaint with index i
 */
double LinIneq::gproj(int i) const {

   return proj[i];

}

/**
 * @return the array containing the projection on the constraints.
 */
double *LinIneq::gproj() {

   return proj;

}

/**
 * @param index the index...
 * @return the value of the constraint labeled with index index
 */
double LinIneq::constraint(int index) const{

   return proj[index] - li[index]->gi()*tr;

}

/**
 * @return the trace of the input TPM scaled with N(N-1)/2
 */
double LinIneq::gtr() const{

   return tr;

}
