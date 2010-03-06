#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * standard constructor\n\n
 * allocates two arrays of dimension n_tp and one of dimension n_ph
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
EIG::EIG(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

}

/**
 * Copy constructor\n\n
 * allocates two arrays of dimension n_tp and one of dimension n_ph, the content
 * of eig_c will be copied into it.
 * @param eig_c input EIG
 */
EIG::EIG(EIG &eig_c){

   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

}
 
/**
 * overload equality operator
 * @param eig_c input EIG that will be copied into this
 */
EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

/**
 * Access to the seperate vectors
 * @param i index 
 * @return corresponding vector
 */
double *EIG::operator[](int i){

   return eig[i];

}

/**
 * Constructor with SUP matrix SZ as input. Allocates two arrays of dimension n_tp and one of dimension n_ph and fills it (this) with the
 * eigevalues of the different SUP blocks. Watch out, the input SUP matrix is destroyed in this operatorion, in it the eigenvectors of the SUP blocks will be stored.
 * @param SZ input SUP
 */
EIG::EIG(SUP &SZ){

   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

   //dit vernietigd de originele matrix!
   for(int i = 0;i < 2;++i)
      (SZ.tpm(i)).diagonalize(eig[i]);

#ifndef PQ

   (SZ.phm()).diagonalize(eig[2]);

#endif

}

/**
 * Destructor, deallocation of the memory.
 */
EIG::~EIG(){

   delete [] eig[0];
   delete [] eig;

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[1][i] << std::endl;

#ifndef PQ

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_ph;++i)
      std::cout << i << "\t" << eig_p.eig[2][i] << std::endl;

#endif

   return output;

}

/**
 */
double EIG::operator()(int block,int index){

   return eig[block][index];

}

/**
 * @return the smallest element in the EIG
 */
double EIG::min(){

   double ward = eig[0][0];

   if(ward > eig[1][0])
      ward = eig[1][0];

#ifndef PQ

   if(ward > eig[2][0])
      ward = eig[2][0];

#endif

   return ward;

}

/**
 * @return the largest element in the EIG
 */
double EIG::max(){

   double ward = eig[0][n_tp - 1];

   if(ward < eig[1][n_tp - 1])
      ward = eig[1][n_tp - 1];

#ifndef PQ

   if(ward < eig[2][n_ph - 1])
      ward = eig[2][n_ph - 1];

#endif

   return ward;

}

/**
 * @return the deviation from the central path as calculed with the logaritmic potential barrier. Should only be used with the eigenvalues as calculated
 * in SUP::center_dev .
 */
double EIG::center_dev(){

   double sum = 0.0;

   for(int i = 0;i < n_tp;++i)
      sum += eig[0][i];

   for(int i = 0;i < n_tp;++i)
      sum += eig[1][i];

   double log_product = 0.0;

   for(int i = 0;i < n_tp;++i)
      log_product += log(eig[0][i]);

   for(int i = 0;i < n_tp;++i)
      log_product += log(eig[1][i]);

#ifndef PQ

   for(int i = 0;i < n_ph;++i)
      sum += eig[2][i];

   for(int i = 0;i < n_ph;++i)
      log_product += log(eig[2][i]);

#endif

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * Returns the deviation from the central path when you take a stepsize alpha in the primal-dual Newton direction (DS,DZ).
 * This is done by calculating the eigenvalues in SUP::line_search, see primal_dual.pdf for more information. \n\n
 * (*this) = eigen_S --> de eigenvalues for the DS step.
 * @param alpha step size along the primal dual direction.
 * @param eigen_Z  --> de eigenvalues for the DZ step.
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter calculated in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter calculated in SUP::line_search
 * @return the deviation from the central path when taking a stepsize alpha in the direction (DS,DZ) from the point (S,Z)
 */
double EIG::centerpot(double alpha,EIG &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig[0][i]);

   return ward;

}
