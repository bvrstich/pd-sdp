#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "SUP.h"
#include "EIG.h"

/**
 * standard constructor\n
 * Alloceert twee TPM matrices en een PHM matrices en zorgt ervoor dat de 
 * daartoe voorziene pointers ernaar verwijzen.
 * @param M aantal sp orbitals
 * @param N aantal deeltjes
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
 * Alloceert twee TPM matrices en een PHM matrices en zorgt ervoor dat de 
 * daartoe voorziene pointers ernaar verwijzen, en kopieert er dan de overeenkomstige matrices uit SZ_c naartoe.
 * @param SZ_c De te kopieren SUP matrix
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
 * Destructor, dealloceerd de voorheen gealloceerde TPM's en PHM.
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
 * @param SZ_pl De SUP matrix die opgeteld moet worden bij this
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
 * @param SZ_pl De SUP matrix die afgetrokken moet worden van this
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
 * @param SZ_c Deze SUP matrix zal gekopieerd worden in this.
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
 * Zet alle getallen in de TPM's en PHM gelijk aan a
 * @param a Het getal
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
 * functie die de TPM blokken teruggeeft
 * @param i = 0 geeft blok 0 terug, = 1 geeft blok 1 terug
 */
TPM &SUP::tpm(int i){

   return *SZ_tp[i];

}

#ifndef PQ

/**
 * functie die de PHM blok teruggeeft, dus eigenlijk blok 2 van de SUP matrix
 */
PHM &SUP::phm(){

   return *SZ_ph;

}

#endif

/**
 * Het initialisatie punt van S, is gewoon u^0: zie primal_dual.pdf voor meer info.
 */
void SUP::init_S(){

   (*SZ_tp[0]).unit();
   (*SZ_tp[1]).Q(1,*SZ_tp[0]);

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
 * Het initialisatie voor Z, zie primal_dual.pdf voor meer info.
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
 * @return aantal deeltjes
 */
int SUP::gN() {

   return N;

}

/**
 * @return aantal sp orbitals
 */
int SUP::gM(){

   return M;

}

/**
 * @return de dimensie de tweedeeltjesruimte
 */
int SUP::gn_tp(){

   return n_tp;

}

#ifndef PQ

/**
 * @return de dimensie van de particle hole ruimte
 */
int SUP::gn_ph(){

   return n_ph;

}

#endif

/**
 * @return de totale dimensie van de SUP matrix
 */
int SUP::gdim(){

   return dim;

}

/**
 * Berekend het inproduct van twee SUP matrices gedefinieerd als Tr(S_1 S_2)
 * @param SZ_i De SUP waarmee het inproduct van this genomen wordt
 * @return double met inproduct in.
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
 * Inverteer de symmetrische, positief definiete SUP matrix.
 * Gebruikt de lapack implementatie van cholesky decompostie dus de matrix MOET
 * positief definiet zijn! Vernietigd de oorspronkelijke matrix.
 */
void SUP::invert(){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->invert();

#ifndef PQ
   
   SZ_ph->invert();

#endif

}

/**
 * Herschaal de SUP matrix met een factor alpha
 * @param alpha Het bewuste getal
 */
void SUP::dscal(double alpha){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->dscal(alpha);

#ifndef PQ
   
   SZ_ph->dscal(alpha);

#endif

}

/**
 * Ortogonale projectie van een algemene SUP matrix diag[ M M_Q M_G] op de U ruimte: diag[M_u Q(M_u) G(M_u)]
 * Voor info kijk naar primal_dual.pdf
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

/**
 * Construeer de D matrix en stop in this, de "metriek" matrix van de hessiaan. Voor meer info zie primal_dual.pdf
 * @param S De primal SUP matrix S
 * @param Z De dual SUP matrix Z
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
 * Neem de positieve of negatieve vierkantswortel uit de SUP matrix. Opgelet, alle matrices worden vernietigd
 * @param option = +1 Neem dan de positieve vierkantswortel, = -1 Neem dan de negatieve vierkantswortel
 */
void SUP::sqrt(int option){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->sqrt(option);

#ifndef PQ

   SZ_ph->sqrt(option);

#endif


}

/**
 * drie maal toepassen van de Matrix::L_map
 */
void SUP::L_map(SUP &map,SUP &object){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->L_map(map.tpm(i),object.tpm(i));

#ifndef PQ

   SZ_ph->L_map(map.phm(),object.phm());

#endif

}

/**
 * Bereken deze SUP plus een constante maal een andere SUP
 * @param alpha Het getal waarmee SZ_p vermenigvuldig wordt
 * @param SZ_p De SUP matrix die alpha maal opgeteld wordt bij this
 */
void SUP::daxpy(double alpha,SUP &SZ_p){

   for(int i = 0;i < 2;++i)
      SZ_tp[i]->daxpy(alpha,SZ_p.tpm(i));

#ifndef PQ
   
   SZ_ph->daxpy(alpha,SZ_p.phm());

#endif

}

/**
 * trace van de SUP matrix, gewoon som van de trace van de TPM 's en PHM (Matrix::trace)
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
 * Ortogonale projectie van een algemene SUP matrix [ M M_Q M_G] op het orthogonaal complement van de U ruimte:
 * Voor info kijk naar primal_dual.pdf
 */
void SUP::proj_C(){

   SUP Z_copy(*this);

   //projecteer op de U ruimte
   Z_copy.proj_U();

   //en het orthogonaal complement nemen:
   *this -= Z_copy;

}

/**
 * algemeen matrixproduct tussen twee SUP matrices, drie maal Matrix::mprod
 * 
 * @param A linkse matrix
 * @param B rechtse matrix
 * @return Het product wordt teruggegeven
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
 * Vul een SUP matrix met TPM matrix: this = diag[tpm  Q(tpm)  G(tpm)]
 * @param tpm De input TPM matrix.
 */
void SUP::fill(TPM &tpm){

   *SZ_tp[0] = tpm;
   SZ_tp[1]->Q(1,tpm);

#ifndef PQ
   
   SZ_ph->G(1,tpm);

#endif

}

/**
 * Implementie van het lineair conjugate gradient algoritme ter oplossing van het duale stelsel\n
 * H(*this) =  B waarin H de hessiaan afbeelding voorstelt.
 * @param B rechterlid van het stelsel
 * @param D SUP matrix die de structuur van de hessiaan afbeelding bepaald. (inverse van die van het primale stelsel)
 * @return return het aantal iteraties dat nodig was om de gewenste nauwkeurigheid te bereiken
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
 * Duale hessiaan afbeelding:\n
 * HB = DBD (dus SUP::L_map) maar dan geprojecteerd op de C ruimte.
 * @param B SUP matrix waarop de hessiaan inwerkt en waarvan de afbeelding wordt opgeslagen in this
 * @param D SUP matrix die de structuur van de hessiaan afbeelding bepaald.
 */
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

/**
 * De U trace van een SUP Z is gedefinieerd als Tr (Zu^0)
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
 * Diagonalisatie van alle blokmatrices d.m.v. Matrix::diagonalize functie, eigenwaarden worden in het EIG object eig gestopt,
 * eigenvectoren zitten zoals bij de Matrix functie in de kolommen van vernietigde oorspronkelijke matrix.
 * @param eig EIG object waar alle eigenwaarden in opgeslaan zitten.
 */
void SUP::diagonalize(EIG &eig){

   SZ_tp[0]->diagonalize(eig[0]);
   SZ_tp[1]->diagonalize(eig[1]);

#ifndef PQ

   SZ_ph->diagonalize(eig[2]);

#endif

}

/**
 * Afwijking van het center via de center_dev potentiaal (logaritmische potentiaal, zie notes) gemeten 
 * Maat voor de afwijking van SZ met de eenheidsmatrix die het moet zijn op het centraal pad\n
 * (*this) = S = primale matrix van het probleem>
 * @param Z = duale matrix van het probleem
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
 * Line search functie die kijkt, voor een gegeven Newton stap (DS,DZ) in de predictor richting,
 * hoe groot de stap "a" is die je kan nemen in die richting voordat de afwijking van het centrum groter is dan max_dev.\n
 * (*this) = DS --> stap van het primale systeem
 * @param DZ stap van het duale systeem
 * @param S Huidige primale SUP
 * @param Z Huidige duale SUP
 * @param max_dev De maximale afwijking van het centrale pad dat je toelaat
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
