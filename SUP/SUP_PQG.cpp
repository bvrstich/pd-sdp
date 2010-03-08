#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "include.h"

/**
 * standard constructor\n\n
 * Calls the SUP_PQ constructor that allocates two TPM matrices, and allocates one extra PHM matrix.
 * @param M number of sp orbitals
 * @param N number particles
 */
SUP_PQG::SUP_PQG(int M,int N) : SUP_PQ(M,N){

   this->n_ph = M*M;
   
   SZ_ph = new PHM(M,N);

   dim += n_ph;

}

/**
 * Copy constructor.\n\n
 * Allocates two TPM matrices and one PHM matrix and copies the content of SZ_c into it.
 * @param SZ_c input SUP_PQG
 */
SUP_PQG::SUP_PQG(SUP_PQG &SZ_c) : SUP_PQ(SZ_c){

   this->n_ph = M*M;

   dim += n_ph;
   
   SZ_ph = new PHM(M,N);

   *SZ_ph = *SZ_c.SZ_ph;

}

/**
 * Destructor, deallocates the memory.
 */
SUP_PQG::~SUP_PQG(){

   delete SZ_ph;

}

/**
 * Overload += operator
 * @param SZ_pl The SUP_PQG matrix that will be added to this
 */
SUP_PQG &SUP_PQG::operator+=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator+=(SZ_pl);

   (*SZ_ph) += (*SZ_pl.SZ_ph);

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP_PQG matrix that will be deducted from this
 */
SUP_PQG &SUP_PQG::operator-=(SUP_PQG &SZ_pl){

   this->SUP_PQ::operator-=(SZ_pl);

   (*SZ_ph) -= (*SZ_pl.SZ_ph);

   return *this;

}

/**
 * overload equality operator
 * @param SZ_c The input SUP_PQG that will be copied into this
 */
SUP_PQG &SUP_PQG::operator=(SUP_PQG &SZ_c){

   this->SUP_PQ::operator=(SZ_pl);

   (*SZ_ph) = (*SZ_c.SZ_ph);

   return *this;

}

/**
 * Make all the numbers in the SUP_PQG equal to a
 * @param a the number
 */
SUP_PQG &SUP_PQG::operator=(double &a){

   this->SUP_PQ::operator=(a);

   (*SZ_ph) = a;

   return *this;

}

/**
 * acces to the PHM block
 * @return The PHM block.
 */
PHM &SUP_PQG::phm(){

   return *SZ_ph;

}

/**
 * Initialisation of the SUP_PQG matix S, is just u^0. (see primal_dual.pdf).
 */
void SUP_PQG::init_S(){

   this->SUP_PQ::init_S();

   SZ_ph->G(1,*SZ_tp[0]);

}

ostream &operator<<(ostream &output,SUP_PQG &SZ_p){

   for(int i = 0;i < 2;++i)
      output << SZ_p.tpm(i) << std::endl;

   output << (*SZ_p.SZ_ph) << std::endl;

   return output;

}

/**
 * Initialisation of the SUP_PQG Z. Must be initialised onto a matrix for which:\n\n
 * Tr (Z u^i) = Tr (H f^i)\n
 * Z >= 0 (positive semidefinite)\n\n
 * see primal_dual.pdf for more information
 */
void SUP_PQG::init_Z(double alpha,TPM &ham,SUP_PQG &u_0){

   (*SZ_tp[0]) = ham;

#ifndef PQ

   *SZ_tp[0] /= 3.0;

#else

   *SZ_tp[0] /= 2.0;

#endif

   SZ_tp[1]->Q(-1,*SZ_tp[0]);

#ifndef PQ
   
   SZ_ph->G(-1,*SZ_tp[0]);

#endif

   //nog een eenheidsmatrix maal constante bijtellen zodat Z positief definiet is:
   this->daxpy(alpha,u_0); 

}

/**
 * @return dimension of ph space
 */
int SUP_PQG::gn_ph(){

   return n_ph;

}

/**
 * Calculates the inproduct of two SUP_PQG matrices defined as Tr (S_1 S_2)
 * @param SZ_i The input SUP_PQG. 
 * @return double with inproduct of this and SZ_i.
 */
double SUP_PQG::ddot(SUP_PQG &SZ_i){

   double ward = this->SUP_PQ::ddot(SZ_i);

   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

   return ward;

}

/**
 * Invert symmetrical positive definite matrix, cholesky decomposition is used so the matrix
 * has to be positive definite!\n
 * original SUP_PQG is destroyed.
 */
void SUP_PQG::invert(){

   this->SUP_PQ::invert();

   SZ_ph->invert();

}

/**
 * Scale the SUP_PQG (this) with factor alpha
 * @param alpha the factor
 */
void SUP_PQG::dscal(double alpha){

   this->SUP_PQ::dscal(alpha);

   SZ_ph->dscal(alpha);

}

/**
 * Take the positive or negative square root of the SUP_PQG, watch out, original SUP_PQG is destroyed (this)
 * @param option = +1 Take the positive square root, = -1 Take the negative square root.
 */
void SUP_PQG::sqrt(int option){

   this->SUP_PQ::sqrt(option);

   SZ_ph->sqrt(option);

}

/**
 * Apply Matrix::L_map onto the different elements in SUP_PQG.
 */
void SUP_PQG::L_map(SUP_PQG &map,SUP_PQG &object){

   this->SUP_PQ::L_map(map,object);

   SZ_ph->L_map(map.phm(),object.phm());

}

/**
 * Calculates this = this + SZ_p *alpha.
 * @param alpha the constant that multiplies SZ_p.
 * @param SZ_p The SUP_PQG matrix that will be added to this.
 */
void SUP_PQG::daxpy(double alpha,SUP_PQG &SZ_p){

   this->SUP_PQ::daxpy(alpha,SZ_p);

   SZ_ph->daxpy(alpha,SZ_p.phm());

}

/**
 * Trace of the SUP_PQG matrix, defined as trace of the seperate blocks.
 */
double SUP_PQG::trace(){

   double ward = this->SUP_PQ::trace();
   
   ward += SZ_ph->trace();

   return ward;

}

/**
 * general matrix product, three times Matrix::mprod. The result is put in this.
 * 
 * @param A left side matrix
 * @param B right side matrix
 */
SUP_PQG &SUP_PQG::mprod(SUP_PQG &A,SUP_PQG &B){

   this->SUP_PQ::mprod(A,B);

   SZ_ph->mprod(A.phm(),B.phm());

   return *this;

}

/**
 * Fill a SUP_PQG matrix with input TPM tpm : this = diag[tpm  Q(tpm)  G(tpm)]
 * @param tpm The input TPM matrix.
 */
void SUP_PQG::fill(TPM &tpm){

   this->SUP_PQ::fill();

   SZ_ph->G(1,tpm);

}

/**
 * Project this onto the space for which Tr (Z u^0) = 0.
 */
void SUP_PQG::proj_U_Tr(){

   //eerst berekenen van de voorfactor voor u^0
   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = this->U_trace();

   double g = (M - N)/(N - 1.0);

   ward /= ( n_tp*(1.0 + q*q) + 
         
         n_ph * (1.0 + g*g) + 2.0*M*g );

   //dan deze factor aftrekken met u^0
   SZ_ph->min_gunit(ward);

   //dan deze factor aftrekken met u^0
   SZ_tp[0]->min_unit(ward);
   SZ_tp[1]->min_qunit(ward);

}

/**
 * @return the U_trace, which is defined as Tr (Z u^0)
 */
double SUP_PQG::U_trace(){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = SZ_tp[0]->trace();

   ward += q*SZ_tp[1]->trace();

   double g = (M - N)/(N - 1.0);

   //skew trace is \sum_{ab} G_{aa;bb}
   ward += g*SZ_ph->trace() + SZ_ph->skew_trace();

   return ward;

}

/**
 * Diagonalisation of all the blockmatrices in SUP_PQG trough Matrix::diagonalize. Eigenvalues are saved in the input EIG_PQG object eig. 
 * Eigenvectors are saved in the original SUP_PQG matrix, so the orignal SUP_PQG this is destroyed.
 * @param eig input EIG_PQG that will containt the eigenvalues after the operation.
 */
EIG_PQG SUP_PQG::diagonalize(){

   EIG_PQG eig(M,N);

   this->SUP_PQ::diagonalize(eig);

   SZ_ph->diagonalize(eig[2]);

   return eig;

}

EIG_PQG *SUP_PQG::get_EIG(){

}
