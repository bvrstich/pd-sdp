#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor:
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param nr the nr of contraints
 * @param option There are predefined constraint, for now there is only option == 0, which means the traceconstraint Tr(Gamma) = N(N-1)/2.
 */
Lineq::Lineq(int M,int N,int nr,int option){

   this->nr = nr;

   this->M = M;
   this->N = N;

   E = new TPM * [nr];
   e = new double [nr];

   for(int i = 0;i < nr;++i)
      E[i] = new TPM(M,N);

   if(option == 0)
      norm_only();
   else{//fill randomly

      for(int i = 0;i < nr;++i){

         E[i]->fill_Random();
         e[i] = (double) rand()/RAND_MAX;

      }

   }

   //orthogonalize the constraints:
   E_ortho = new TPM * [nr];
   e_ortho = new double [nr];

   for(int i = 0;i < nr;++i)
      E_ortho[i] = new TPM(M,N);

   orthogonalize();//speaks for itself, doesn't it?

   //now construct the u_0
   u_0 = new SUP(M,N);

   u_0->tpm(0) = 0;

   for(int i = 0;i < nr;++i)
      (u_0->tpm(0)).daxpy(e_ortho[i],*E_ortho[i]);

   u_0->fill();

   u_0_norm = u_0->ddot(*u_0);

}

/**
 * copy constructor:
 * @param lineq the Lineq object that will be copied into this
 */
Lineq::Lineq(const Lineq &lineq){

   this->nr = lineq.gnr();

   this->M = lineq.gM();
   this->N = lineq.gN();

   E = new TPM * [nr];
   e = new double [nr];

   E_ortho = new TPM * [nr];
   e_ortho = new double [nr];

   //copy the that which needs to be copied
   for(int i = 0;i < nr;++i){

      E[i] = new TPM(lineq.gE(i)); 
      e[i] = lineq.ge(i);

      //no need to orthogonalize, just copy!
      E_ortho[i] = new TPM(lineq.gE_ortho(i)); 
      e_ortho[i] = lineq.ge_ortho(i);

   }

   //now construct the u_0
   u_0 = new SUP(lineq.gu_0());

   //and its norm.
   u_0_norm = lineq.gu_0_norm();

}

/**
 * destructor
 */
Lineq::~Lineq(){

   for(int i = 0;i < nr;++i){

      delete E[i];
      delete E_ortho[i];

   }

   delete [] E;
   delete [] E_ortho;

   delete [] e;
   delete [] e_ortho;

   delete u_0;

}

/**
 * @return nr of particles
 */
int Lineq::gN() const {

   return N;

}

/**
 * @return nr of constraints
 */
int Lineq::gnr() const {

   return nr;

}

/**
 * @return nr of sp orbitals
 */
int Lineq::gM() const {

   return M;

}

/**
 * access to the individual constraint TPM's
 * @param i the index
 * @return The E TPM on index i.
 */
TPM &Lineq::gE(int i) const {

   return *E[i];

}

/**
 * access to the individual constraint values
 * @param i the index
 * @return the e values on index i: e[i] or something
 */
double &Lineq::ge(int i) const{

   return e[i];
   
}

/**
 * access to the individual orthogonalized constraint TPM's: private function
 * @param i the index
 * @return The E_ortho TPM on index i.
 */
TPM &Lineq::gE_ortho(int i) const {

   return *E_ortho[i];

}

/**
 * access to the individual orthogonalized constraint values: private function
 * @param i the index
 * @return the e values on index i: e_ortho[i] or something
 */
double &Lineq::ge_ortho(int i) const{

   return e_ortho[i];
   
}

ostream &operator<<(ostream &output,Lineq &lineq_p){

   output << "first print the constraint matrices:";
   output << endl;
   output << endl;

   for(int i = 0;i < lineq_p.gnr();++i){

      output << "constraint nr :" << i << endl;
      output << endl;

      output << lineq_p.gE(i);

   }

   output << endl;
   output << endl;
   output << "the desired values are:" << endl;
   output << endl;

   for(int i = 0;i < lineq_p.gnr();++i)
      output << i << "\t" << lineq_p.ge(i) << endl;

   return output;

}

/**
 * orthogonalize the constraints, will take E and e, and construct E_ortho and e_ortho with them.
 */
void Lineq::orthogonalize(){

   //construct the overlapmatrix of the E's
   Matrix S(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         S(i,j) = E[i]->ddot(*E[j]);

   S.symmetrize();

   //take the inverse square root
   S.sqrt(-1);

   //make the orthogonal ones:
   for(int i = 0;i < nr;++i){

      *E_ortho[i] = 0;
      e_ortho[i] = 0;

      for(int j = 0;j < nr;++j){

         E_ortho[i]->daxpy(S(i,j),*E[j]);
         e_ortho[i] += S(i,j) * e[j];

      }

   }

}

/**
 * access to the u_0
 * @return u_0
 */
SUP &Lineq::gu_0() const{

   return *u_0;

}

/**
 * @return the norm of u_0
 */
double Lineq::gu_0_norm() const{

   return u_0_norm;

}

/**
 * fills the constraint matrices with the norm constraint only
 */
void Lineq::norm_only(){

   *E[0] = 0;

   for(int i = 0;i < E[0]->gn();++i)
      (*E[0])(i,i) = 1.0;

   e[0] = N*(N - 1.0)/2.0;

}
