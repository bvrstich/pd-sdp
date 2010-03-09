#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * standard constructor, allocates the memory for the eigenvalues of a SUP object
 * @param M dimension of sp space
 * @param N nr of particles
 */
EIG::EIG(int M,int N){

   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

   this->dim = 2*n_tp;

#ifdef __G_CON
   
   this->n_ph = M*M;

   dim += n_ph;

#endif
   
   eig = new double [dim];

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(EIG &eig_c){

   this->N = eig_c.N;
   this->M = eig_c.M;

   this->n_tp = M*(M - 1)/2;

   this->dim = 2*n_tp;

#ifdef __G_CON
   
   this->n_ph = M*M;

   dim += n_ph;

#endif
   
   eig = new double [dim];

   int inc = 1;

   dcopy_(&dim,eig_c.eig,&inc,eig,&inc);

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig,&inc,eig,&inc);

   return *this;

}

/** 
 * get pointer to the eigenvalues
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned, i == 2 eigenvalues of the G block, etc.
 * @return array of eigenvalues
 */
double *EIG::operator[](int i){

   return eig + i*n_tp;

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param SZ input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &SZ){

   //first allocate the memory
   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

   this->dim = 2*n_tp;

#ifdef __G_CON
   
   this->n_ph = M*M;

   dim += n_ph;

#endif
   
   eig = new double [dim];

   //then diagonalize the SUP
   for(int i = 0;i < 2;++i)
      (SZ.tpm(i)).diagonalize(eig + i*n_tp);

#ifdef __G_CON

   (SZ.phm()).diagonalize(eig + 2*n_tp);

#endif

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete [] eig;

}

ostream &operator<<(ostream &output,EIG &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p(0,i) << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p(1,i) << std::endl;

#ifdef __G_CON

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_ph;++i)
      std::cout << i << "\t" << eig_p(2,i) << std::endl;

#endif

   return output;

}

/**
 * acces to the numbers
 * @param block == 0, get element "index" from P block, == 1 get element index from Q block, == 2 get element index from G block 
 * @param index which element in block you want
 * @return eig[block][index]
 */
double EIG::operator()(int block,int index){

   return eig[block*n_tp + index];

}

/**
 * @return nr of particles
 */
int EIG::gN(){

   return N;

}

/**
 * @return dimension of sp space
 */
int EIG::gM(){

   return M;

}

/**
 * @return dimension of tp space
 */
int EIG::gn_tp(){

   return n_tp;

}

#ifdef __G_CON

/**
 * @return dimension of tp space
 */
int EIG::gn_ph(){

   return n_ph;

}

#endif

/**
 * @return total dimension of the EIG object
 */
int EIG::gdim(){

   return dim;

}


/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min(){

   //lowest eigenvalue of P block
   double ward = eig[0];

   //lowest eigenvalue of Q block
   if(ward > eig[n_tp])
      ward = eig[n_tp];

#ifdef __G_CON

   //lowest eigenvalue of G block
   if(ward > eig[2*n_tp])
      ward = eig[2*n_tp];

#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max(){

   //maximum of P block
   double ward = eig[n_tp - 1];

   //maximum of Q block
   if(ward < eig[2*n_tp - 1])
      ward = eig[2*n_tp - 1];

#ifdef __G_CON

   //maximum of G block
   if(ward < eig[2*n_tp + n_ph - 1])
      ward = eig[2*n_tp + n_ph - 1];

#endif

   return ward;

}

/**
 * @return The deviation of the central path as calculated with the logarithmic barrierfunction, the EIG object is calculated
 * in SUP::center_dev.
 */
double EIG::center_dev(){

   double sum = 0.0;

   for(int i = 0;i < dim;++i)
      sum += eig[i];

   double log_product = 0.0;

   for(int i = 0;i < dim;++i)
      log_product += log(eig[i]);

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * @return the deviation of the central path measured trough the logarithmic potential barrier (see primal_dual.pdf), when you take a stepsize alpha from
 * the point (S,Z) in the primal dual newton direction (DS,DZ), for which you have calculated the generalized eigenvalues eigen_S and eigen_Z in SUP::line_search.
 * (*this) = eigen_S --> generalized eignevalues for the DS step
 * @param alpha the stepsize
 * @param eigen_Z --> generalized eignevalues for the DS step
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter calculated in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter calculated in SUP::line_search
 */
double EIG::centerpot(double alpha,EIG &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eig[i]);

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig[i]);

   return ward;

}
