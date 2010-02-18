#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "EIG.h"
#include "SUP.h"
#include "lapack.h"

//constructor:
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

//copy constructor:
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
 
EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

double *EIG::operator[](int i){

   return eig[i];

}

//constructor met initialisatie door SUP matrix:
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

//destructor
EIG::~EIG(){

   delete [] eig[0];
   delete [] eig;

}

//friend function! output stream operator overloaded
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

double EIG::operator()(int block,int index){

   return eig[block][index];

}

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

double EIG::centerpot(double alpha,EIG &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig[0][i]);

   return ward;

}
