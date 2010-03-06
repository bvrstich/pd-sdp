#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "include.h"

/**
 * standard constructor\n\n
 * Allocates two TPM matrices and one PHM matrix.
 * @param M number of sp orbitals
 * @param N number particles
 */
SUP::SUP(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   this->dim = 2*n_tp;

#ifndef PQ

   this->n_ph = M*M;
   
   SZ_ph = new PHM(M,N);

   dim += n_ph;

#endif

}

/**
 * Copy constructor.\n\n
 * Allocates two TPM matrices and one PHM matrix and copies the content of SZ_c into it.
 * @param SZ_c input SUP
 */
SUP::SUP(SUP &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;
   this->dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifndef PQ

   this->n_ph = M*M;

   dim += n_ph;
   
   SZ_ph = new PHM(M,N);

   *SZ_ph = *SZ_c.SZ_ph;

#endif

}

/**
 * Destructor, deallocates the memory.
 */
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

#ifndef PQ
   
   delete SZ_ph;

#endif

}

/**
 * Overload += operator
 * @param SZ_pl The SUP matrix that will be added to this
 */
SUP &SUP::operator+=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

#ifndef PQ
   
   (*SZ_ph) += (*SZ_pl.SZ_ph);

#endif

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP matrix that will be deducted from this
 */
SUP &SUP::operator-=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

#ifndef PQ
   
   (*SZ_ph) -= (*SZ_pl.SZ_ph);

#endif

   return *this;

}

/**
 * overload equality operator
 * @param SZ_c The input SUP that will be copied into this
 */
SUP &SUP::operator=(SUP &SZ_c){

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifndef PQ
   
   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

   return *this;

}

/**
 * Make all the numbers in the SUP equal to a
 * @param a the number
 */
SUP &SUP::operator=(double &a){

   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

#ifndef PQ
   
   (*SZ_ph) = a;

#endif

   return *this;

}

/**
 * acces to the seperate TPM blocks
 * @param i = 0 returns block 0 (the \Gamma block), = 1 returns block 1 (the Q block)
 * @return The TPM element corresponding to the index i
 * 
 */
TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

#ifndef PQ

/**
 * acces to the PHM block
 * @return The PHM block.
 */
PHM &SUP::phm(){

   return *SZ_ph;

}

#endif

/**
 * Initialisation of the SUP matix S, is just u^0. (see primal_dual.pdf).
 */
void SUP::init_S(){

   SZ_tp[0]->unit();
   SZ_tp[1]->Q(1,*SZ_tp[0]);

#ifndef PQ
   
   SZ_ph->G(1,*SZ_tp[0]);

#endif

}

ostream &operator<<(ostream &output,SUP &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

#ifndef PQ
   
   output << std::endl;
   output << (*SZ_p.SZ_ph);

#endif

   return output;

}

/**
 * Initialisation of the SUP Z. Must be initialised onto a matrix for which:\n\n
 * Tr (Z u^i) = Tr (H f^i)\n
 * Z >= 0 (positive semidefinite)\n\n
 * see primal_dual.pdf for more information
 */
void SUP::init_Z(double alpha,TPM &ham,SUP &u_0){

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
 * @return number of particles
 */
int SUP::gN() {

   return N;

}

/**
 * @return number of sp orbitals
 */
int SUP::gM(){

   return M;

}

/**
 * @return dimension of tp space
 */
int SUP::gn_tp(){

   return n_tp;

}

#ifndef PQ

/**
 * @return dimension of ph space
 */
int SUP::gn_ph(){

   return n_ph;

}

#endif

/**
 * @return the total dimension of the SUP matrix.
 */
int SUP::gdim(){

   return dim;

}

/**
 * Calculates the inproduct of two SUP matrices defined as Tr (S_1 S_2)
 * @param SZ_i The input SUP. 
 * @return double with inproduct of this and SZ_i.
 */
double SUP::ddot(SUP &SZ_i){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

#ifndef PQ
   
   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

#endif

   return ward;

}

/**
 * Invert symmetrical positive definite matrix, cholesky decomposition is used so the matrix
 * has to be positive definite!\n
 * original SUP is destroyed.
 */
void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifndef PQ
   
   SZ_ph->invert();

#endif

}

/**
 * Scale the SUP (this) with factor alpha
 * @param alpha the factor
 */
void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifndef PQ
   
   SZ_ph->dscal(alpha);

#endif

}

/**
 * Orthogonal projection of a general SUP matrix diag[ M M_Q M_G] onto the U-space: diag[M_u Q(M_u) G(M_u)].\n
 * see primal_dual.pdf for more information.
 */
void SUP::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) + G(M_G) in O stoppen
   TPM O(*SZ_tp[0]);

   SZ_tp[0]->Q(1,*SZ_tp[1]);

   O += *SZ_tp[0];

#ifndef PQ
   
   SZ_tp[0]->G(1,*SZ_ph);

   O += *SZ_tp[0];

#endif

   //dan de inverse overlapmatrix hierop laten inwerken en in this[0] stoppen
   SZ_tp[0]->S(-1,O);

   //en de Q hiervan in this[1]
   SZ_tp[1]->Q(1,*SZ_tp[0]);

#ifndef PQ

   //en de G hiervan in in phm van SZ
   SZ_ph->G(1,*SZ_tp[0]);

#endif
   
   //Nu is de projectie op de u^\alpha's gebeurd.
   //Nu nog de projectie op de u^i's: dus component langs u^0 eruit halen
   this->proj_U_Tr();

}

/**
 * Construct the D matrix and put it in this, D is the metric matrix of the hessian, see primal_dual.pdf.
 * @param S The primal SUP matrix S
 * @param Z The dual SUP matrix Z
 */
void SUP::D(SUP &S,SUP &Z){

   //positieve vierkantswortel uit Z
   SUP Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP hulp(M,N);

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Take the positive or negative square root of the SUP, watch out, original SUP is destroyed (this)
 * @param option = +1 Take the positive square root, = -1 Take the negative square root.
 */
void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifndef PQ

   SZ_ph->sqrt(option);

#endif


}

/**
 * Apply Matrix::L_map onto the different elements in SUP.
 */
void SUP::L_map(SUP &map,SUP &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifndef PQ

   SZ_ph->L_map(map.phm(),object.phm());

#endif

}

/**
 * Calculates this = this + SZ_p *alpha.
 * @param alpha the constant that multiplies SZ_p.
 * @param SZ_p The SUP matrix that will be added to this.
 */
void SUP::daxpy(double alpha,SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifndef PQ
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

}

/**
 * Trace of the SUP matrix, defined as trace of the seperate blocks.
 */
double SUP::trace(){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

#ifndef PQ
   
   ward += SZ_ph->trace();

#endif

   return ward;

}

/**
 * Ortogonal projection of a general SUP matrix [ M M_Q M_G] onto the ortogonal complement of the U-space. (see primal_dual.pdf)
 */
void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

/**
 * general matrix product, three times Matrix::mprod. The result is put in this.
 * 
 * @param A left side matrix
 * @param B right side matrix
 */
SUP &SUP::mprod(SUP &A,SUP &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

#ifndef PQ

   SZ_ph->mprod(A.phm(),B.phm());

#endif

   return *this;

}

/**
 * Fill a SUP matrix with input TPM tpm : this = diag[tpm  Q(tpm)  G(tpm)]
 * @param tpm The input TPM matrix.
 */
void SUP::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifndef PQ
   
   SZ_ph->G(1,tpm);

#endif

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton equation:
 * H(*this) =  B in which H is the Hessian map SUP::H
 * @param B right hand side of the equation
 * @param D SUP matrix that determines the structure of the hessian map ( as calculated in SUP::D ) 
 * @return the number of iterations needed to converge.
 */
int SUP::solve(SUP &B,SUP &D){

   SUP HB(M,N);
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP r(B);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:

   int cg_iter = 0;

   while(rr > 1.0e-3){

      ++cg_iter;

      HB.H(B,D);

      ward = rr/B.ddot(HB);

      //delta += ward*b
      this->daxpy(ward,B);

      //r -= ward*HB
      r.daxpy(-ward,HB);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      B.dscal(rr/rr_old);

      //dan r er bijtellen
      B += r;

   }
   
   return cg_iter;

}

/**
 * Duale hessian map:\n\n
 * HB = DBD (dus SUP::L_map) projected onto the C-space (as calculated with SUP::proj_C )
 * @param B SUP matrix onto which the Hessian map works, image is put in this.
 * @param D SUP matrix that defines the structure of the hessian map. (metric)
 */
void SUP::H(SUP &B,SUP &D){

   this->L_map(D,B);

   this->proj_C();

}

/**
 * Project this onto the space for which Tr (Z u^0) = 0.
 */
void SUP::proj_U_Tr(){

   //eerst berekenen van de voorfactor voor u^0
   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = this->U_trace();

#ifdef PQ

   ward /= ( n_tp*(1.0 + q*q) );

#else

   double g = (M - N)/(N - 1.0);

   ward /= ( n_tp*(1.0 + q*q) + 
         
         n_ph * (1.0 + g*g) + 2.0*M*g );

   //dan deze factor aftrekken met u^0
   SZ_ph->min_gunit(ward);

#endif
   
   //dan deze factor aftrekken met u^0
   SZ_tp[0]->min_unit(ward);
   SZ_tp[1]->min_qunit(ward);

}

/**
 * @return the U_trace, which is defined as Tr (Z u^0)
 */
double SUP::U_trace(){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = SZ_tp[0]->trace();

   ward += q*SZ_tp[1]->trace();

#ifndef PQ

   double g = (M - N)/(N - 1.0);

   //skew trace is \sum_{ab} G_{aa;bb}
   ward += g*SZ_ph->trace() + SZ_ph->skew_trace();

#endif

   return ward;

}

/**
 * Diagonalisation of all the blockmatrices in SUP trough Matrix::diagonalize. Eigenvalues are saved in the input EIG object eig. 
 * Eigenvectors are saved in the original SUP matrix, so the orignal SUP this is destroyed.
 * @param eig input EIG that will containt the eigenvalues after the operation.
 */
void SUP::diagonalize(EIG &eig){

   SZ_tp[0]->diagonalize(eig[0]);
   SZ_tp[1]->diagonalize(eig[1]);

#ifndef PQ

   SZ_ph->diagonalize(eig[2]);

#endif

}

/**
 * Deviation from the central path.\n\n
 * (*this) = S = primal SUP matrix of the problem
 * @param Z = dual SUP matrix of the problem
 * @return Measure of the deviation from the central path through the logaritmic potential function (See primal_dua.pdf).
 */
double SUP::center_dev(SUP &Z){

   SUP sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   EIG eig(SZ);

   return eig.center_dev();

}

/**
 * Line search function that checks, for a given primal dual Newton direction (DS,DZ) in the predictor direction,
 * how large a step in this direction you can take before deviating more then max_dev from the central path.
 * (*this) = DS --> primal Newton direction.
 * @param DZ dual Newton direction.
 * @param S Current primal SUP
 * @param Z Current dual SUP
 * @param max_dev The maximal deviation from the central path you want to allow.
 */
double SUP::line_search(SUP &DZ,SUP &S,SUP &Z,double max_dev){

   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP hulp(M,N);
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG eigen_S(hulp);

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG eigen_Z(hulp);

   //nog c_S en c_Z uitrekenen:
   double pd_gap = S.ddot(Z);

   //c_S = Tr (DS Z)/Tr (SZ)
   double c_S = this->ddot(Z)/pd_gap;

   //c_Z = Tr (S DZ)/Tr (SZ)
   double c_Z = S.ddot(DZ)/pd_gap;

   //waar zitten de singulariteiten: tot waar mag ik zoeken?
   double a_max = -1.0/eigen_S.min();
   double b_max = -1.0/eigen_Z.min();

   //a_max is de waarde tot waar ik zal zoeken:
   if(b_max < a_max)
      a_max = b_max;

   double a = 0.0;
   double b = a_max;

   double c = (a + b)/2.0;

   //bissectiemethode om stapgrootte te bepalen:
   while(b - a > 1.0e-5){

      c = (a + b)/2.0;

      if( (center_dev + eigen_S.centerpot(c,eigen_Z,c_S,c_Z) - max_dev) < 0.0 )
         a = c;
      else
         b = c;

   }

   return c;

}
