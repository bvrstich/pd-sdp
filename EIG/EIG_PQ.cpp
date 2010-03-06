#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "EIG/EIG_PQ.h"
#include "SUP.h"
#include "lapack.h"

/**
 * standard constructor\n\n
 * allocates two arrays of dimension n_tp and one of dimension n_ph
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
EIG_PQ::EIG_PQ(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

}

/**
 * Copy constructor\n\n
 * allocates two arrays of dimension n_tp and one of dimension n_ph, the content
 * of eig_c will be copied into it.
 * @param eig_c input EIG_PQ
 */
EIG_PQ::EIG_PQ(EIG_PQ &eig_c){

   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

}
 
/**
 * overload equality operator
 * @param eig_c input EIG_PQ that will be copied into this
 */
EIG_PQ &EIG_PQ::operator=(EIG_PQ &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

/** 
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned.
 * @return array of eigenvalues
 */
double *EIG_PQ::operator[](int i){

   return eig[i];

}

/**
 * Constructor with SUP matrix SZ as input. Allocates two arrays of dimension n_tp and one of dimension n_ph and fills it (this) with the
 * eigevalues of the different SUP blocks. Watch out, the input SUP matrix is destroyed in this operatorion, in it the eigenvectors of the SUP blocks will be stored.
 * @param SZ input SUP
 */
EIG_PQ::EIG_PQ(SUP_PQ &SZ){

   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

   //dit vernietigd de originele matrix!
   for(int i = 0;i < 2;++i)
      (SZ.tpm(i)).diagonalize(eig[i]);

}

/**
 * Destructor, deallocation of the memory.
 */
EIG_PQ::~EIG_PQ(){

   delete [] eig[0];
   delete [] eig;

}

ostream &operator<<(ostream &output,const EIG_PQ &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[1][i] << std::endl;

   return output;

}

/**
 * acces to the numbers
 * @param block == 0, get element "index" from P block, == 1 get element index from Q block 
 * @param index which element in block you want
 * @return eig[block][index]
 */
double EIG_PQ::operator()(int block,int index){

   return eig[block][index];

}

/**
 * @return the smallest element in the EIG_PQ
 */
double EIG_PQ::min(){

   double ward = eig[0][0];

   if(ward > eig[1][0])
      ward = eig[1][0];

   return ward;

}

/**
 * @return the largest element in the EIG_PQ
 */
double EIG_PQ::max(){

   double ward = eig[0][n_tp - 1];

   if(ward < eig[1][n_tp - 1])
      ward = eig[1][n_tp - 1];

   return ward;

}

/**
 * @return the deviation from the central path as calculed with the logaritmic potential barrier. Should only be used with the eigenvalues as calculated
 * in SUP::center_dev .
 */
double EIG_PQ::center_dev(){

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
double EIG_PQ::centerpot(double alpha,EIG_PQ &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig[0][i]);

   return ward;

}
