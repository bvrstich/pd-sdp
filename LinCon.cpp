#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * Constructor of a LinCon object
 * @param M The constraint matrix
 * @param N the minimal projection
 */
LinCon::LinCon(int M,int N){

   I_c = new TPM(M,N);

   I_c_bar = new SPM(M,N);

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param lc_copy The LinCon object to be copied
 */
LinCon::LinCon(const LinCon &lc_copy){

   I_c = new TPM(lc_copy.gI());

   I_c_bar = new SPM(lc_copy.gI_bar());

   i_c = lc_copy.gi();

   this->M = lc_copy.gM();
   this->N = lc_copy.gN();

}

/**
 * destructor
 */
LinCon::~LinCon(){

   delete I_c;

   delete I_c_bar;

}

/**
 * @return the constraint TPM object
 */
TPM &LinCon::gI() const{

   return *I_c;

}

/**
 * @return the partially trace constraint, the SPM object I_c_bar.
 */
SPM &LinCon::gI_bar() const{

   return *I_c_bar;

}

/**
 * @return The minimal projection
 */
double LinCon::gi() const{

   return i_c;

}

/**
 * set the constraint value
 * @param i the value that the minimal projection will be set to.
 */
void LinCon::si(double i){

   i_c = i;

}

/**
 * set the constraint Matrix, warning, first set the value i_c, because otherwise the shift will be wrong
 * @param I the input constraint Matrix
 */
void LinCon::sI(const TPM &I){

   *I_c = I;

   for(int i = 0;i < I.gn();++i)//shift it for convenience
      (*I_c)(i,i) -= i_c/(N*(N - 1.0));

   I_c_bar->bar(I);

}

ostream &operator<<(ostream &output,const LinCon &lc_p){

   cout << "The minimal projection:\t" << lc_p.gi() << endl;
   cout << endl;

   cout << "The shifted Constraint matrix:" << endl;
   cout << endl;

   cout << lc_p.gI() << endl;

   return output;

}

/**
 * @return nr of sp orbitals
 */
int LinCon::gM() const{

   return M;

}

/**
 * @return nr of particles
 */
int LinCon::gN() const{

   return N;

}

/**
 * construct a diagonal T constraint for diagonal element index
 * @param index the index of the diagonal element for which the constraint matrix will be constructed
 */
void LinCon::diag_T(int index){

   PPHM lincon(M,N);

   lincon = 0.0;
   lincon(index,index) = 1.0;

   I_c->T(lincon);

   I_c_bar->bar(*I_c);

   i_c = 0.0;

}

/**
 * construct the spin matrix as the spin matrix
 * @param spin the value of the spinconstraint
 */
void LinCon::spincon(double spin){

   I_c->set_S_2();

   for(int i = 0;i < I_c->gn();++i)
      (*I_c)(i,i) -= spin/(N*(N - 1.0));

   I_c_bar->bar(*I_c);

   i_c = spin;

}

/**
 * fills the constraints object randomly
 */
void LinCon::fill_Random(){

   i_c = rand()/RAND_MAX;

   I_c->fill_Random();

   for(int i = 0;i < I_c->gn();++i)
      (*I_c)(i,i) -= i_c/(N*(N - 1.0));

   I_c_bar->bar(*I_c);

}
