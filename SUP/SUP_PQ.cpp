#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "SUP/SUP_PQ.h"
#include "EIG/EIG_PQ.h"

/**
 * standard constructor\n\n
 * Allocates the two TPM matrices .
 * @param M number of sp orbitals
 * @param N number particles
 */
SUP_PQ::SUP_PQ(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   this->dim = 2*n_tp;

}

/**
 * Copy constructor.\n\n
 * Allocates two the TPM matrices copies the content of SZ_c into it.
 * @param SZ_c input SUP_PQ
 */
SUP_PQ::SUP_PQ(SUP_PQ &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;
   this->dim = 2*n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

}

/**
 * Destructor, deallocates the memory.
 */
SUP_PQ::~SUP_PQ(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

}

/**
 * Overload += operator
 * @param SZ_pl The SUP_PQ matrix that will be added to this
 */
SUP_PQ &SUP_PQ::operator+=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * Overload -= operator
 * @param SZ_pl The SUP_PQ matrix that will be deducted from this
 */
SUP_PQ &SUP_PQ::operator-=(SUP_PQ &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

   return *this;

}

/**
 * overload equality operator
 * @param SZ_c The input SUP_PQ that will be copied into this
 */
SUP_PQ &SUP_PQ::operator=(SUP_PQ &SZ_c){

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

   return *this;

}

/**
 * Make all the numbers in the SUP_PQ equal to the number a
 * @param a the number
 */
SUP_PQ &SUP_PQ::operator=(double &a){

   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

   return *this;

}

/**
 * acces to the seperate TPM blocks
 * @param i = 0 returns block 0 (the \Gamma block), = 1 returns block 1 (the Q block)
 * @return The TPM element corresponding to the index i
 * 
 */
TPM &SUP_PQ::tpm(int i){

   return *SZ_tp[i];

}

/**
 * Initialisation of the SUP_PQ matix S, is just u^0. (see primal_dual.pdf).
 */
void SUP_PQ::init_S(){

   SZ_tp[0]->unit();
   SZ_tp[1]->Q(1,*SZ_tp[0]);

}

ostream &operator<<(ostream &output,SUP_PQ &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

   return output;

}

/**
 * Initialisation of the SUP_PQ Z. Must be initialised onto a matrix for which:\n\n
 * Tr (Z u^i) = Tr (H f^i)\n
 * Z >= 0 (positive semidefinite)\n\n
 * see primal_dual.pdf for more information
 */
void SUP_PQ::init_Z(double alpha,TPM &ham,SUP_PQ &u_0){

   (*SZ_tp[0]) = ham;

   *SZ_tp[0] /= 2.0;

   SZ_tp[1]->Q(-1,*SZ_tp[0]);

   //nog een eenheidsmatrix maal constante bijtellen zodat Z positief definiet is:
   this->daxpy(alpha,u_0); 

}

/**
 * @return number of particles
 */
int SUP_PQ::gN() {

   return N;

}

/**
 * @return number of sp orbitals
 */
int SUP_PQ::gM(){

   return M;

}

/**
 * @return dimension of tp space
 */
int SUP_PQ::gn_tp(){

   return n_tp;

}

/**
 * @return the total dimension of the SUP_PQ matrix.
 */
int SUP_PQ::gdim(){

   return dim;

}

/**
 * Calculates the inproduct of two SUP_PQ matrices defined as Tr (S_1 S_2)
 * @param SZ_i The input SUP_PQ. 
 * @return double with inproduct of this and SZ_i.
 */
double SUP_PQ::ddot(SUP_PQ &SZ_i){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

   return ward;

}

/**
 * Invert symmetrical positive definite matrix, cholesky decomposition is used so the matrix
 * has to be positive definite!\n
 * original SUP_PQ is destroyed.
 */
void SUP_PQ::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

}

/**
 * Scale the SUP_PQ (this) with factor alpha
 * @param alpha the factor
 */
void SUP_PQ::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

}

/**
 * Orthogonal projection of a general SUP_PQ matrix diag[ M M_Q M_G] onto the U-space: diag[M_u Q(M_u) G(M_u)].\n
 * see primal_dual.pdf for more information.
 */
void SUP_PQ::proj_U(){
  
   //eerst M_Gamma + Q(M_Q) in O stoppen
   TPM O(M,N);

   O.collaps(*this);

   //dan de inverse overlapmatrix hierop laten inwerken en in hulp stoppen
   TPM hulp(M,N);

   hulp.S(-1,O);

   //en this vullen met deze TPM
   this->fill(hulp);

   //Nu is de projectie op de u^\alpha's gebeurd.
   //Nu nog de projectie op de u^i's: dus component langs u^0 eruit halen
   this->proj_U_Tr();

}

/**
 * Construct the D matrix and put it in this, D is the metric matrix of the hessian, see primal_dual.pdf.
 * @param S The primal SUP_PQ matrix S
 * @param Z The dual SUP_PQ matrix Z
 */
void SUP_PQ::D(SUP_PQ &S,SUP_PQ &Z){

   //positieve vierkantswortel uit Z
   SUP_PQ Z_copy(Z);

   Z_copy.sqrt(1);

   //links en rechts vermenigvuldigen met wortel Z
   SUP_PQ hulp(M,N);

   hulp.L_map(Z_copy,S);

   hulp.sqrt(1);

   //negatieve vierkantswortel uit Z
   Z_copy = Z;

   Z_copy.sqrt(-1);

   //en links en rechts hulp vermenigvuldigen met die wortel, en in this steken:
   this->L_map(Z_copy,hulp);

}

/**
 * Take the positive or negative square root of the SUP_PQ, watch out, original SUP_PQ is destroyed (this)
 * @param option = +1 Take the positive square root, = -1 Take the negative square root.
 */
void SUP_PQ::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

}

/**
 * Apply Matrix::L_map onto the different elements in SUP_PQ.
 */
void SUP_PQ::L_map(SUP_PQ &map,SUP_PQ &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

}

/**
 * Calculates this = this + SZ_p *alpha.
 * @param alpha the constant that multiplies SZ_p.
 * @param SZ_p The SUP_PQ matrix that will be added to this.
 */
void SUP_PQ::daxpy(double alpha,SUP_PQ &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

}

/**
 * Trace of the SUP_PQ matrix, defined as trace of the seperate blocks.
 */
double SUP_PQ::trace(){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

   return ward;

}

/**
 * Ortogonal projection of a general SUP_PQ matrix [ M M_Q ] onto the ortogonal complement of the U-space. (see primal_dual.pdf)
 */
void SUP_PQ::proj_C(){

   SUP_PQ Z_copy(*this);

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
SUP_PQ &SUP_PQ::mprod(SUP_PQ &A,SUP_PQ &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

   return *this;

}

/**
 * Fill a SUP_PQ matrix with input TPM tpm : this = diag[tpm  Q(tpm)]
 * @param tpm The input TPM matrix.
 */
void SUP_PQ::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

}

/**
 * Implementation of the linear conjugate gradient algorithm for the solution of the dual Newton equation:
 * H(*this) =  B in which H is the Hessian map SUP_PQ::H
 * @param B right hand side of the equation
 * @param D SUP_PQ matrix that determines the structure of the hessian map ( as calculated in SUP_PQ::D ) 
 * @return the number of iterations needed to converge.
 */
int SUP_PQ::solve(SUP_PQ &B,SUP_PQ &D){

   SUP_PQ HB(M,N);
   HB.H(*this,D);

   B -= HB;

   //de r initialiseren op B - H DZ
   SUP_PQ r(B);

   double rr = r.ddot(r);
   double rr_old,ward;

   int cg_iter = 0;

   while(rr > 1.0e-5){

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
 * Dual hessian map:\n\n
 * HB = DBD (dus SUP_PQ::L_map) projected onto the C-space (as calculated with SUP_PQ::proj_C )
 * @param B SUP_PQ matrix onto which the Hessian map works, image is put in this.
 * @param D SUP_PQ matrix that defines the structure of the hessian map. (metric)
 */
void SUP_PQ::H(SUP_PQ &B,SUP_PQ &D){

   this->L_map(D,B);

   this->proj_C();

}

/**
 * Project (*this) onto the space for which Tr (Z u^0) = 0.
 */
void SUP_PQ::proj_U_Tr(){

   //eerst berekenen van de voorfactor voor u^0
   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = this->U_trace();

   ward /= ( n_tp*(1.0 + q*q) );

   //dan deze factor aftrekken met u^0
   SZ_tp[0]->min_unit(ward);
   SZ_tp[1]->min_qunit(ward);

}

/**
 * @return the U_trace, which is defined as Tr (Z u^0)
 */
double SUP_PQ::U_trace(){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   double ward = SZ_tp[0]->trace();

   ward += q*SZ_tp[1]->trace();

   return ward;

}

/**
 * Diagonalisation of all the blockmatrices in SUP_PQ trough Matrix::diagonalize. Eigenvalues are saved in the EIG_PQ object eig. 
 * Eigenvectors are saved in the original SUP_PQ matrix (*this), so the original SUP_PQ this is destroyed.
 * @return EIG_PQ that will contain the eigenvalues of this after the operation.
 */
EIG_PQ SUP_PQ::diagonalize(){

   EIG_PQ eig(M,N);

   SZ_tp[0]->diagonalize(eig[0]);
   SZ_tp[1]->diagonalize(eig[1]);

   return eig;

}

/**
 * Diagonalisation of all the blockmatrices in SUP_PQ trough Matrix::diagonalize. Eigenvalues are saved in a EIG_PQ object eig. 
 * Eigenvectors are saved in the original SUP_PQ matrix (*this), so the original SUP_PQ this is destroyed.\n
 * Watch out: a EIG_PQ object here is dynamically allocated and has to be deleted explicitally after calling this function.
 * @return pointer to the EIG_PQ that will contain the eigenvalues of this after the operation.
 */
EIG_PQ *SUP_PQ::get_EIG(){

   EIG_PQ *eig = new EIG_PQ(*this);

   return eig;

}

/**
 * Deviation from the central path.\n\n
 * (*this) = S = primal SUP_PQ matrix of the problem
 * @param Z = dual SUP_PQ matrix of the problem
 * @return Measure of the deviation from the central path through the logaritmic potential function (See primal_dua.pdf).
 */
double SUP_PQ::center_dev(SUP_PQ &Z){

   SUP_PQ sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP_PQ SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   return (SZ.diagonalize()).center_dev();


}

/**
 * Line search function that checks, for a given primal dual Newton direction (DS,DZ) in the predictor direction,
 * how large a step in this direction you can take before deviating more then max_dev from the central path.
 * (*this) = DS --> primal Newton direction.
 * @param DZ dual Newton direction.
 * @param S Current primal SUP_PQ
 * @param Z Current dual SUP_PQ
 * @param max_dev The maximal deviation from the central path you want to allow.
 */
double SUP_PQ::line_search(SUP_PQ &DZ,SUP_PQ &S,SUP_PQ &Z,double max_dev){

   //eerst de huidige deviatie van het centraal pad nemen:
   double center_dev = S.center_dev(Z);

   //eigenwaarden zoeken van S^{-1/2} DS S^{-1/2} en Z^{-1/2} DZ Z^{-1/2}

   //kopieer S in de zogeheten wortel:
   SUP_PQ wortel(S);

   //maak negatieve vierkantswortel uit S
   wortel.sqrt(-1);

   //de L_map
   SUP_PQ hulp(M,N);
   hulp.L_map(wortel,*this);

   //eigenwaarden in eigen_S stoppen
   EIG_PQ *eigen_S = hulp.get_EIG();

   //nu idem voor Z
   wortel = Z;

   wortel.sqrt(-1);

   hulp.L_map(wortel,DZ);

   EIG_PQ *eigen_Z = hulp.get_EIG();

   //nog c_S en c_Z uitrekenen:
   double pd_gap = S.ddot(Z);

   //c_S = Tr (DS Z)/Tr (SZ)
   double c_S = this->ddot(Z)/pd_gap;

   //c_Z = Tr (S DZ)/Tr (SZ)
   double c_Z = S.ddot(DZ)/pd_gap;

   //waar zitten de singulariteiten: tot waar mag ik zoeken?
   double a_max = -1.0/eigen_S->min();
   double b_max = -1.0/eigen_Z->min();

   //a_max is de waarde tot waar ik zal zoeken:
   if(b_max < a_max)
      a_max = b_max;

   double a = 0.0;
   double b = a_max;

   double c = (a + b)/2.0;

   //bissectiemethode om stapgrootte te bepalen:
   while(b - a > 1.0e-5){

      c = (a + b)/2.0;

      if( (center_dev + eigen_S->centerpot(c,*eigen_Z,c_S,c_Z) - max_dev) < 0.0 )
         a = c;
      else
         b = c;

   }

   //nog allocated memory deleten
   delete eigen_S;
   delete eigen_Z;

   return c;

}
