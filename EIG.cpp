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

   eig_tp = new double [2*n_tp];

   this->dim = 2*n_tp;

#ifdef __G_CON
   
   this->n_ph = M*M;

   eig_ph = new double [n_ph];

   dim += n_ph;

#endif

#ifdef __T1_CON
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   eig_dp = new double [n_dp];

   dim += n_dp;

#endif

#ifdef __T2_CON

   this->n_pph = M*M*(M - 1)/2;

   eig_pph = new double [n_pph];

   dim += n_pph;

#endif

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

   eig_tp = new double [2*n_tp];

   this->dim = 2*n_tp;

   int inc = 1;

   dcopy_(&dim,eig_c.eig_tp,&inc,eig_tp,&inc);

#ifdef __G_CON
   
   this->n_ph = M*M;

   eig_ph = new double [n_ph];

   dim += n_ph;

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

#endif
   
#ifdef __T1_CON
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   eig_dp = new double [n_dp];

   dim += n_dp;

   dcopy_(&n_dp,eig_c.eig_dp,&inc,eig_dp,&inc);

#endif

#ifdef __T2_CON
   
   this->n_pph = M*M*(M - 1)/2;

   eig_pph = new double [n_pph];

   dim += n_pph;

   dcopy_(&n_pph,eig_c.eig_pph,&inc,eig_pph,&inc);

#endif


}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   int n = 2*n_tp;

   dcopy_(&n,eig_c.eig_tp,&inc,eig_tp,&inc);

#ifdef __G_CON

   dcopy_(&n_ph,eig_c.eig_ph,&inc,eig_ph,&inc);

#endif

#ifdef __T1_CON

   dcopy_(&n_dp,eig_c.eig_dp,&inc,eig_dp,&inc);

#endif

#ifdef __T2_CON

   dcopy_(&n_pph,eig_c.eig_pph,&inc,eig_pph,&inc);

#endif

   return *this;

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

   eig_tp = new double [2*n_tp];

   this->dim = 2*n_tp;

   //diagonalize
   for(int i = 0;i < 2;++i)
      (SZ.tpm(i)).diagonalize(eig_tp + i*n_tp);

#ifdef __G_CON
   
   this->n_ph = M*M;

   eig_ph = new double [n_ph];

   dim += n_ph;

   //diagonalize
   (SZ.phm()).diagonalize(eig_ph);

#endif

#ifdef __T1_CON
   
   this->n_dp = M*(M - 1)*(M - 2)/6;

   eig_dp = new double [n_dp];

   dim += n_dp;

   //diagonalize
   (SZ.dpm()).diagonalize(eig_dp);

#endif

#ifdef __T2_CON

   this->n_pph = M*M*(M - 1)/2;

   eig_pph = new double [n_pph];

   dim += n_pph;

   //diagonalize
   (SZ.pphm()).diagonalize(eig_pph);
 
#endif

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete [] eig_tp;

#ifdef __G_CON
   
   delete [] eig_ph;

#endif

#ifdef __T1_CON

   delete [] eig_dp;

#endif

#ifdef __T2_CON

   delete [] eig_pph;

#endif

}

ostream &operator<<(ostream &output,EIG &eig_p){

   for(int i = 0;i < eig_p.gn_tp();++i)
      std::cout << i << "\t" << eig_p.tpm(0)[i] << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.gn_tp();++i)
      std::cout << i << "\t" << eig_p.tpm(1)[i] << std::endl;

#ifdef __G_CON

   std::cout << std::endl;

   for(int i = 0;i < eig_p.gn_ph();++i)
      std::cout << i << "\t" << eig_p.phm()[i] << std::endl;

#endif

#ifdef __T1_CON

   std::cout << std::endl;

   for(int i = 0;i < eig_p.gn_dp();++i)
      std::cout << i << "\t" << eig_p.dpm()[i] << std::endl;

#endif

#ifdef __T2_CON

   std::cout << std::endl;

   for(int i = 0;i < eig_p.gn_pph();++i)
      std::cout << i << "\t" << eig_p.pphm()[i] << std::endl;

#endif

   return output;

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

/** 
 * get pointer to the eigenvalues of the TPM blocks P and Q
 * @param i == 0, the eigenvalues of the P block will be returned, i == 1, the eigenvalues of the Q block will be returned
 * @return array of eigenvalues
 */
double *EIG::tpm(int i){

   return eig_tp + i*n_tp;

}


#ifdef __G_CON

/**
 * @return dimension of ph space
 */
int EIG::gn_ph(){

   return n_ph;

}

/** 
 * get pointer to the eigenvalues of the PHM block of the SUP
 * @return array of eigenvalues
 */
double *EIG::phm(){

   return eig_ph;

}

#endif

#ifdef __T1_CON

/**
 * @return dimension of dp space
 */
int EIG::gn_dp(){

   return n_dp;

}

/** 
 * get pointer to the eigenvalues of the DPM block of the SUP
 * @return array of eigenvalues
 */
double *EIG::dpm(){

   return eig_dp;

}

#endif

#ifdef __T2_CON

/**
 * @return dimension of pph space
 */
int EIG::gn_pph(){

   return n_pph;

}

/** 
 * get pointer to the eigenvalues of the PPHM block of the SUP
 * @return array of eigenvalues
 */
double *EIG::pphm(){

   return eig_pph;

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
   double ward = eig_tp[0];

   //lowest eigenvalue of Q block
   if(ward > eig_tp[n_tp])
      ward = eig_tp[n_tp];

#ifdef __G_CON

   //lowest eigenvalue of G block
   if(ward > eig_ph[0])
      ward = eig_ph[0];

#endif

#ifdef __T1_CON

   //lowest eigenvalue of T1 block
   if(ward > eig_dp[0])
      ward = eig_dp[0];

#endif

#ifdef __T2_CON

   //lowest eigenvalue of T2 block
   if(ward > eig_pph[0])
      ward = eig_pph[0];

#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max(){

   //maximum of P block
   double ward = eig_tp[n_tp - 1];

   //maximum of Q block
   if(ward < eig_tp[2*n_tp - 1])
      ward = eig_tp[2*n_tp - 1];

#ifdef __G_CON

   //maximum of G block
   if(ward < eig_ph[n_ph - 1])
      ward = eig_ph[n_ph - 1];

#endif

#ifdef __T1_CON

   //maximum of T1 block
   if(ward < eig_dp[n_dp - 1])
      ward = eig_dp[n_dp - 1];

#endif

#ifdef __T2_CON

   //maximum of T2 block
   if(ward < eig_pph[n_pph - 1])
      ward = eig_pph[n_pph - 1];

#endif

   return ward;

}

/**
 * @return The deviation of the central path as calculated with the logarithmic barrierfunction, the EIG object is calculated
 * in SUP::center_dev.
 */
double EIG::center_dev(){

   double sum = 0.0;

   for(int i = 0;i < 2*n_tp;++i)
      sum += eig_tp[i];

   double log_product = 0.0;

   for(int i = 0;i < 2*n_tp;++i)
      log_product += log(eig_tp[i]);

#ifdef __G_CON

   for(int i = 0;i < n_ph;++i)
      sum += eig_ph[i];

   for(int i = 0;i < n_ph;++i)
      log_product += log(eig_ph[i]);

#endif

#ifdef __T1_CON

   for(int i = 0;i < n_dp;++i)
      sum += eig_dp[i];

   for(int i = 0;i < n_dp;++i)
      log_product += log(eig_dp[i]);

#endif

#ifdef __T2_CON

   for(int i = 0;i < n_pph;++i)
      sum += eig_pph[i];

   for(int i = 0;i < n_pph;++i)
      log_product += log(eig_pph[i]);

#endif

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * @return the deviation of the central path measured trough the logarithmic potential barrier (see primal_dual.pdf), when you take a stepsize alpha from
 * the point (S,Z) in the primal dual newton direction (DS,DZ), for which you have calculated the generalized eigenvalues eigen_S and eigen_Z in SUP::line_search.
 * (*this) = eigen_S --> generalized eignevalues for the DS step
 * @param alpha the stepsize
 * @param eigen_Z --> generalized eigenvalues for the DS step
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter calculated in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter calculated in SUP::line_search
 */
double EIG::centerpot(double alpha,EIG &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < 2*n_tp;++i)
      ward -= log(1.0 + alpha*eig_tp[i]);

   for(int i = 0;i < 2*n_tp;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig_tp[i]);

#ifdef __G_CON

   for(int i = 0;i < n_ph;++i)
      ward -= log(1.0 + alpha*eig_ph[i]);

   for(int i = 0;i < n_ph;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig_ph[i]);

#endif

#ifdef __T1_CON

   for(int i = 0;i < n_dp;++i)
      ward -= log(1.0 + alpha*eig_dp[i]);

   for(int i = 0;i < n_dp;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig_dp[i]);

#endif

#ifdef __T2_CON

   for(int i = 0;i < n_pph;++i)
      ward -= log(1.0 + alpha*eig_pph[i]);

   for(int i = 0;i < n_pph;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig_pph[i]);

#endif

   return ward;

}
