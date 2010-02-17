#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "SUP.h"
#include "EIG.h"

//constructor
SUP::SUP(int M,int N){

   this->M = M;
   this->N = N;
   this->n_tp = M*(M - 1)/2;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

#ifndef PQ

   this->n_ph = M*M;
   
   SZ_ph = new PHM(M,N);

#endif

}

//copy constructor
SUP::SUP(SUP &SZ_c){

   this->M = SZ_c.M;
   this->N = SZ_c.N;
   this->n_tp = SZ_c.n_tp;

   SZ_tp = new TPM * [2];

   for(int i = 0;i < 2;++i)
      SZ_tp[i] = new TPM(M,N);

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifndef PQ

   this->n_ph = M*M;
   
   SZ_ph = new PHM(M,N);

   *SZ_ph = *SZ_c.SZ_ph;

#endif

}

//destructor
SUP::~SUP(){

   for(int i = 0;i < 2;++i)
      delete SZ_tp[i];

   delete [] SZ_tp;

#ifndef PQ
   
   delete SZ_ph;

#endif

}

SUP &SUP::operator+=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) += (*SZ_pl.SZ_tp[i]);

#ifndef PQ
   
   (*SZ_ph) += (*SZ_pl.SZ_ph);

#endif

   return *this;

}

SUP &SUP::operator-=(SUP &SZ_pl){

   for(int i = 0;i < 2;++i)
      (*SZ_tp[i]) -= (*SZ_pl.SZ_tp[i]);

#ifndef PQ
   
   (*SZ_ph) -= (*SZ_pl.SZ_ph);

#endif

   return *this;

}

//overload equality operator
SUP &SUP::operator=(SUP &SZ_c){

   (*SZ_tp[0]) = (*SZ_c.SZ_tp[0]);
   (*SZ_tp[1]) = (*SZ_c.SZ_tp[1]);

#ifndef PQ
   
   (*SZ_ph) = (*SZ_c.SZ_ph);

#endif

   return *this;

}

SUP &SUP::operator=(double &a){

   (*SZ_tp[0]) = a;
   (*SZ_tp[1]) = a;

#ifndef PQ
   
   (*SZ_ph) = a;

#endif

   return *this;

}

TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

#ifndef PQ

PHM &SUP::phm(){

   return *SZ_ph;

}

#endif

void SUP::init_S(){

   (*SZ_tp[0]).unit();
   (*SZ_tp[1]).Q(1,*SZ_tp[0]);

#ifndef PQ
   
   SZ_ph->G(1,*SZ_tp[0]);

#endif

}

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,SUP &SZ_p){

   output << (*SZ_p.SZ_tp[0]) << std::endl;
   output << (*SZ_p.SZ_tp[1]);

#ifndef PQ
   
   output << std::endl;
   output << (*SZ_p.SZ_ph);

#endif

   return output;

}

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

int SUP::gN() {

   return N;

}

int SUP::gM(){

   return M;

}

int SUP::gn_tp(){

   return n_tp;

}

#ifndef PQ

int SUP::gn_ph(){

   return n_ph;

}

#endif

double SUP::ddot(SUP &SZ_i){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->ddot(*SZ_i.SZ_tp[i]);

#ifndef PQ
   
   ward += SZ_ph->ddot(*SZ_i.SZ_ph);

#endif

   return ward;

}

void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifndef PQ
   
   SZ_ph->invert();

#endif

}

void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifndef PQ
   
   SZ_ph->dscal(alpha);

#endif

}

//projecteer SUP matrix *this op de U ruimte
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
   (*SZ_tp[0]).S(-1,O);

   //en de Q hiervan in this[1]
   (*SZ_tp[1]).Q(1,*SZ_tp[0]);

#ifndef PQ

   //en de G hiervan in in phm van SZ
   SZ_ph->G(1,*SZ_tp[0]);

#endif
   
   //Nu is de projectie op de u^\alpha's gebeurd.
   //Nu nog de projectie op de u^i's: dus component langs u^0 eruit halen
   this->proj_U_Tr();

}

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

void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifndef PQ

   SZ_ph->sqrt(option);

#endif


}

void SUP::L_map(SUP &map,SUP &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifndef PQ

   SZ_ph->L_map(map.phm(),object.phm());

#endif

}

void SUP::daxpy(double alpha,SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifndef PQ
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

}

double SUP::trace(){

   double ward = 0.0;

   for(int i = 0;i < 2;++i)
      ward += SZ_tp[i]->trace();

#ifndef PQ
   
   ward += SZ_ph->trace();

#endif

   return ward;

}

void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

SUP &SUP::mprod(SUP &A,SUP &B){

   for(int i= 0;i < 2;++i)
      SZ_tp[i]->mprod(A.tpm(i),B.tpm(i));

#ifndef PQ

   SZ_ph->mprod(A.phm(),B.phm());

#endif

   return *this;

}

void SUP::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifndef PQ
   
   SZ_ph->G(1,tpm);

#endif

}

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

void SUP::H(SUP &B,SUP &D){

   this->L_map(D,B);

   this->proj_C();

}

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

void SUP::diagonalize(EIG &eig){

   SZ_tp[0]->diagonalize(eig[0]);
   SZ_tp[1]->diagonalize(eig[1]);

#ifndef PQ

   SZ_ph->diagonalize(eig[2]);

#endif

}

//afwijking van het center via de potentiaal gemeten
//this = S, afwijking van center in combinatie met Z
double SUP::center_dev(SUP &Z){

   SUP sqrt_S(*this);

   sqrt_S.sqrt(1);

   SUP SZ(M,N);
   SZ.L_map(sqrt_S,Z);

   EIG eig(SZ);

   return eig.center_dev();

}
