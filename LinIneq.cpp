#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int LinIneq::nr;

LinCon **LinIneq::li;

/**
 * Static function that initialializes the static members: it allocates the LinCon objects and 
 * initializes them.
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param nr_in nr of contraint matrices
 */
void LinIneq::init(int M,int N,int nr_in){

   nr = nr_in;

   //allocate
   li = new LinCon * [nr];

   for(int i = 0;i < nr;++i)
      li[i] = new LinCon(M,N);

   //fill
   for(int i = 0;i < nr;++i)
      li[i]->fill_Random();

}

/**
 * delete the statics
 */
void LinIneq::clean(){

   for(int i = 0;i < nr;++i)
      delete li[i];

   delete [] li;

}

/**
 * constructor of a LinIneq object
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
LinIneq::LinIneq(int M,int N){

   this->M = M;
   this->N = N;

   proj = new double [nr];

}

/**
 * copy constructor
 * @param li_copy The LinIneq object to be copied
 */
LinIneq::LinIneq(const LinIneq &li_copy){

   this->M = li_copy.gM();
   this->N = li_copy.gN();

   this->tr = li_copy.gtr();

   proj = new double [nr];

   for(int i = 0;i < nr;++i)
      proj[i] = li_copy.gproj(i);

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   delete [] proj;

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
