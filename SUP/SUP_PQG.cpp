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
 * Orthogonal projection of a general SUP_PQG matrix diag[ M M_Q M_G] onto the U-space: diag[M_u Q(M_u) G(M_u)].\n
 * see primal_dual.pdf for more information.
 */
void SUP_PQG::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) + G(M_G) in O stoppen (collapsen)
   TPM O(M,N);

   O.collaps(*this);

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
 * @param S The primal SUP_PQG matrix S
 * @param Z The dual SUP_PQG matrix Z
 */
void SUP_PQG::D(SUP_PQG &S,SUP_PQG &Z){

   //positieve vierkantswortel uit Z
   SUP_PQG Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP_PQG hulp(M,N);

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Take the positive or negative square root of the SUP_PQG, watch out, original SUP_PQG is destroyed (this)
 * @param option = +1 Take the positive square root, = -1 Take the negative square root.
 */
void SUP_PQG::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifndef PQ

   SZ_ph->sqrt(option);

#endif


}

/**
 * Apply Matrix::L_map onto the different elements in SUP_PQG.
 */
void SUP_PQG::L_map(SUP_PQG &map,SUP_PQG &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifndef PQ

   SZ_ph->L_map(map.phm(),object.phm());

#endif

}

/**
 * Calculates this = this + SZ_p *alpha.
 * @param alpha the constant that multiplies SZ_p.
 * @param SZ_p The SUP_PQG matrix that will be added to this.
 */
void SUP_PQG::daxpy(double alpha,SUP_PQG &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifndef PQ
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

}

/**
 * Trace of the SUP_PQG matrix, defined as trace of the seperate blocks.
 */
double SUP_PQG::trace(){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

#ifndef PQ
   
   ward += SZ_ph->trace();

#endif

   return ward;

}

/**
 * Ortogonal projection of a general SUP_PQG matrix [ M M_Q M_G] onto the ortogonal complement of the U-space. (see primal_dual.pdf)
 */
void SUP_PQG::proj_C(){

   SUP_PQG Z_copy(*this);

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
SUP_PQG &SUP_PQG::mprod(SUP_PQG &A,SUP_PQG &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

#ifndef PQ

   SZ_ph->mprod(A.phm(),B.phm());

#endif

   return *this;

}

/**
 * Fill a SUP_PQG matrix with input TPM tpm : this = diag[tpm  Q(tpm)  G(tpm)]
 * @param tpm The input TPM matrix.
 */
void SUP_PQG::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifndef PQ
   
   SZ_ph->G(1,tpm);

#endif

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton equation:
 * H(*this) =  B in which H is the Hessian map SUP_PQG::H
 * @param B right hand side of the equation
 * @param D SUP_PQG matrix that determines the structure of the hessian map ( as calculated in SUP_PQG::D ) 
 * @return the number of iterations needed to converge.
 */
int SUP_PQG::solve(SUP_PQG &B,SUP_PQG &D){

   SUP_PQG HB(M,N);
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP_PQG r(B);

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
 * HB = DBD (dus SUP_PQG::L_map) projected onto the C-space (as calculated with SUP_PQG::proj_C )
 * @param B SUP_PQG matrix onto which the Hessian map works, image is put in this.
 * @param D SUP_PQG matrix that defines the structure of the hessian map. (metric)
 */
void SUP_PQG::H(SUP_PQG &B,SUP_PQG &D){

   this->L_map(D,B);

   this->proj_C();

}

/**
 * Project this onto the space for which Tr (Z u^0) = 0.
 */
void SUP_PQG::proj_U_Tr(){

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
double SUP_PQG::U_trace(){

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
 * Diagonalisation of all the blockmatrices in SUP_PQG trough Matrix::diagonalize. Eigenvalues are saved in the input EIG_PQG object eig. 
 * Eigenvectors are saved in the original SUP_PQG matrix, so the orignal SUP_PQG this is destroyed.
 * @param eig input EIG_PQG that will containt the eigenvalues after the operation.
 */
void SUP_PQG::diagonalize(EIG_PQG &eig){

   SZ_tp[0]->diagonalize(eig[0]);
   SZ_tp[1]->diagonalize(eig[1]);

#ifndef PQ

   SZ_ph->diagonalize(eig[2]);

#endif

}

/**
 * Deviation from the central path.\n\n
 * (*this) = S = primal SUP_PQG matrix of the problem
 * @param Z = dual SUP_PQG matrix of the problem
 * @return Measure of the deviation from the central path through the logaritmic potential function (See primal_dua.pdf).
 */
double SUP_PQG::center_dev(SUP_PQG &Z){

   SUP_PQG sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP_PQG SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   EIG_PQG eig(SZ);

   return eig.center_dev();

}

/**
 * Line search function that checks, for a given primal dual Newton direction (DS,DZ) in the predictor direction,
 * how large a step in this direction you can take before deviating more then max_dev from the central path.
 * (*this) = DS --> primal Newton direction.
 * @param DZ dual Newton direction.
 * @param S Current primal SUP_PQG
 * @param Z Current dual SUP_PQG
 * @param max_dev The maximal deviation from the central path you want to allow.
 */
double SUP_PQG::line_search(SUP_PQG &DZ,SUP_PQG &S,SUP_PQG &Z,double max_dev){

   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP_PQG wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP_PQG hulp(M,N);
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG_PQG eigen_S(hulp);

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG_PQG eigen_Z(hulp);

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
