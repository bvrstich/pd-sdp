#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::endl;

#include "SPM.h"
#include "TPM.h"

#ifndef PQ

#include "PHM.h"

#endif

#include "SUP.h"
#include "lapack.h"

#include "DPM.h"

int TPM::counter = 0;

int **TPM::t2s;
int **TPM::s2t;

/**
 * standard constructor, maakt Matrix aan met dimensie M*(M - 1)/2
 * als counter == 0 dan alloceerd de constructor het geheugen voor lijsten t2s en s2t en initialiseerd deze:
 * @param M aantal sp orbitals
 * @param N aantal deeltjes
 */
TPM::TPM(int M,int N) : Matrix(M*(M - 1)/2) {

   this->N = N;
   this->M = M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

/**
 * copy constructor, maakt Matrix aan met dimensie M*(M - 1)/2 en kopieerd er tpm_c in
 * als counter == 0 dan alloceerd de constructor het geheugen voor lijsten t2s en s2t en initialiseerd deze:
 * @param tpm_c matrix die gekopieerd moet worden
 */
TPM::TPM(TPM &tpm_c) : Matrix(tpm_c){

   this->N = tpm_c.N;
   this->M = tpm_c.M;
   this->n = M*(M - 1)/2;

   if(counter == 0){

      //allocatie van sp2tp
      s2t = new int * [M];
      s2t[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2t[i] = s2t[i - 1] + M;

      //allocatie van tp2sp
      t2s = new int * [n];

      for(int i = 0;i < n;++i)
         t2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j){

            s2t[i][j] = teller;

            t2s[teller][0] = i;
            t2s[teller][1] = j;

            ++teller;

         }

      for(int i = 0;i < M;++i)
         for(int j = i + 1;j < M;++j)
            s2t[j][i] = s2t[i][j];

   }

   ++counter;

}

/**
 * Destructor: als counter == 1 dan dealloceerd de destructor het geheugen voor lijsten t2s en s2t.
 */
TPM::~TPM(){

   if(counter == 1){

      delete [] s2t[0];
      delete [] s2t;

      for(int i = 0;i < n;++i)
         delete [] t2s[i];

      delete [] t2s;

   }

   --counter;

}

/**
 * toegang tot de getallen in de matrix door gebruik te maken van de sp indices,
 * er wordt rekenging gehouden met de antisymmetrie.\n
 * vb. TPM(a,b,c,d) = - TPM(b,a,c,d) = -TPM(a,b,d,c) = TPM(b,a,d,c)
 * @param a eerste sp index, vormt samen met b de eerste tp index
 * @param b tweede sp index, vormt samen met a de eerste tp index
 * @param c derde sp index, vormt samen met d de tweede tp index
 * @param d vierde sp index, vormt samen met c de tweede tp index
 */
double TPM::operator()(int a,int b,int c,int d) const{

   if( (a == b) || (c == d) )
      return 0;
   else{

      int i = s2t[a][b];
      int j = s2t[c][d];

      int phase = 1;

      if(a > b)
         phase *= -1;
      if(c > d)
         phase *= -1;

      return phase*(*this)(i,j);

   }

}

ostream &operator<<(ostream &output,TPM &tpm_p){

   for(int i = 0;i < tpm_p.n;++i)
      for(int j = 0;j < tpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return aantal deeltjes
 */
int TPM::gN(){

   return N;

}

/**
 * @return aantal sp orbitals
 */
int TPM::gM(){

   return M;

}

/**
 * @return de dimensie van de matrix en van de tweedeeltjesruimte
 */
int TPM::gn(){

   return n;

}

/**
 * Maakt de gereduceerde tweedeeltjes hubbard hamiltoniaan aan en stop hem in this
 * @param U de sterkte van de on site repulsie (U > 0) of attractie (U < 0)
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp orbitals

   double ward = 1.0/(N - 1.0);

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0;

         //eerst hopping
         if( (a == c) && ( ( (b + 2)%M == d ) || ( b == (d + 2)%M ) ) )
            (*this)(i,j) -= ward;

         if( (b == c) && ( ( (a + 2)%M == d ) || ( a == (d + 2)%M ) ) )
            (*this)(i,j) += ward;

         if( (b == d) && ( ( (a + 2)%M == c ) || ( a == (c + 2)%M ) ) )
            (*this)(i,j) -= ward;

         //on site interaction
         if( (a % 2) == 0 && (c % 2) == 0 )
            if(a == (b - 1) && c == (d - 1) && a == c)
               (*this)(i,j) += U;

      }

   }

   this->symmetrize();

}

/**
 * De Q afbeelding
 * @param option = 1, gewone Q afbeelding , = -1 inverse Q afbeelding
 * @param tpm_d De TPM waarvan de Q-like afbeelding genomen wordt en in this gestoken wordt
 */
void TPM::Q(int option,TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * De Q-like afbeelding: zie primal-dual.pdf voor meer info (vorm: Q(A,B,C)(TPM) )
 * @param option = 1, gewone Q afbeelding , = -1 inverse Q afbeelding
 * @param A voorfactor van two particle stuk afbeelding
 * @param B voorfactor van no particle stuk afbeelding
 * @param C voorfactor van single particle stuk afbeelding
 * @param tpm_d De TPM waarvan de Q-like afbeelding genomen wordt en in this gestoken wordt
 */
void TPM::Q(int option,double A,double B,double C,TPM &tpm_d){

   if(option == -1){

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm(M,N);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   //construct de spm met schaling C
   spm.constr(C,tpm_d);

   for(int i = 0;i < n;++i){

      int a = t2s[i][0];
      int b = t2s[i][1];

      for(int j = i;j < n;++j){

         int c = t2s[j][0];
         int d = t2s[j][1];

         (*this)(i,j) = A*tpm_d(i,j);

         if(i == j)
            (*this)(i,i) += ward;

         if(a == c)
            (*this)(i,j) -= spm(b,d);

         if(b == c)
            (*this)(i,j) += spm(a,d);

         if(b == d)
            (*this)(i,j) -= spm(a,c);

      }
   }

   this->symmetrize();

}
/**
 * initialiseer this op de eenheidsmatrix met trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = N*(N - 1.0)/(2.0*n);

   for(int i = 0;i < n;++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < n;++j)
         (*this)(i,j) = (*this)(j,i) = 0.0;

   }

}

/**
 * orthogonale projectie op traceless ruimte
 */
void TPM::proj_Tr(){

   double ward = (this->trace())/(double)n;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= ward;

}

/**
 * Primale hessiaan afbeelding:\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3\n
 * met D_1,D_2 en D_3 de P,Q en G blokken van de SUP D. (Indien gecompileerd wordt met optie PQ wordt het G stuk weggelaten)
 * @param b TPM matrix waarop de hessiaan inwerkt en waarvan de afbeelding wordt opgeslagen in this
 * @param D SUP matrix die de structuur van de hessiaan afbeelding bepaald.
 */
void TPM::H(TPM &b,SUP &D){

   this->L_map(D.tpm(0),b);

   //maak Q(b)
   TPM Qb(M,N);
   Qb.Q(1,b);

   TPM hulp(M,N);

   hulp.L_map(D.tpm(1),Qb);

   Qb.Q(1,hulp);

   *this += Qb;

#ifndef PQ

   //maak G(b)
   PHM Gb(M,N);
   Gb.G(1,b);

   PHM hulpje(M,N);

   hulpje.L_map(D.phm(),Gb);

   hulp.G(1,hulpje);

   *this += hulp;

#endif

   this->proj_Tr();

}

/**
 * Implementie van het lineair conjugate gradient algoritme ter oplossing van het primale stelsel\n
 * H(*this) =  b waarin H de hessiaan afbeelding voorstelt.
 * @param b rechterlid van het stelsel
 * @param D SUP matrix die de structuur van de hessiaan afbeelding bepaald.
 * @return return het aantal iteraties dat nodig was om de gewenste nauwkeurigheid te bereiken
 */
int TPM::solve(TPM &b,SUP &D){

   *this = 0;

   //de r initialiseren op b
   TPM r(b);

   double rr = r.ddot(r);
   double rr_old,ward;

   //nog de Hb aanmaken ook, maar niet initialiseren:
   TPM Hb(M,N);

   int cg_iter = 0;

   while(rr > 1.0e-5){

      ++cg_iter;

      Hb.H(b,D);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe r_norm maken
      rr_old = rr;
      rr = r.ddot(r);

      //eerst herschalen van b
      b.dscal(rr/rr_old);

      //dan r er bijtellen
      b += r;

   }

   return cg_iter;

}


#ifndef PQ
/**
 * De G afbeelding die een PHM object afbeeld op een TPM object.
 * @param option = 1 dan wordt G_down uitgevoerd, = -1 dan wordt G^{-1}_up uitgevoerd
 * @param phm input PHM die afgebeeld wordt op this
 */
void TPM::G(int option,PHM &phm){

   SPM spm(M,N);

   if(option == 1)
      spm.constr(1.0/(N - 1.0),phm);
   else
      spm.constr(1.0/(M - N + 1.0),phm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = phm(b,d,c,a) - phm(a,d,c,b) - phm(b,c,d,a) + phm(a,c,d,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

         if(b == c)
            (*this)(i,j) -= spm(a,d);

         if(a == c)
            (*this)(i,j) += spm(b,d);

      }

   }

   //nog schalen met 4 voor G^{-1}
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

#endif

/**
 * Overlapmatrix afbeelding, is eigenlijk een Q-like afbeelding waarvoor ik de 
 * parameters a,b en c heb berekend in primal-dual.pdf. Aangezien het een Q-like afbeelding is hebben we dus onmiddelijk ook de inverse
 * overlapmatrix-afbeelding.
 * @param option = 1 directe overlapmatrixafbeelding , = -1 inverse overlapmatrix afbeelding
 * @param tpm_d de input TPM die afgebeeld wordt op this
 */
void TPM::S(int option,TPM &tpm_d){

#ifdef PQ

   double a = 2.0;
   double b = (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   double c = (2.0*N - M)/((N - 1.0)*(N - 1.0));

#else

   double a = 6.0;
   double b = (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   double c = (4.0*N - 2.0*M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

   this->Q(option,a,b,c,tpm_d);

}

/**
 * Trek van this de eenheidsmatrix * een constante af:\n
 * this -= scale* 1
 * @param scale de constante waarmee de eenheidsmarix vermenigvuldigt wordt
 */
void TPM::min_unit(double scale){

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * Trek van this - de Q-afbeelding van de eenheidsmatrix  * een constante - af:\n
 * this -= scale* Q(1)
 * @param scale de constante waarmee de eenheidsmarix vermenigvuldigt wordt
 */
void TPM::min_qunit(double scale){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   scale *= q;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

/**
 * calculate the trace of one pair of sp indices of a DPM an put in (*this):\n\n
 * TPM(a,b,d,e) = \sum_{c} DPM(a,b,c,d,e,c)
 * @param dpm input DPM
 */
void TPM::bar(DPM &dpm){

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = 0;j < n;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += dpm(a,b,l,c,d,l);

      }
   }

   this->symmetrize();

}

/**
 * map a DPM (dpm) on a TPM (*this) with a T1 map, (Q-like map), watch out for the inverse
 * up map, when M = 2*N it is singular! So don't use it!:
 * @param option = +1 T1_down , =-1 inverse T1_up
 * @param dpm The input DPM
 */
void TPM::T(int option,DPM &dpm){

   TPM tpm(M,N);
   tpm.bar(dpm);

   if(option == 1){

      double a = 1;
      double b = 1.0/(3.0*N*(N - 1.0));
      double c = 0.5/(N - 1.0);

      this->Q(1,a,b,c,tpm);

   }
   else{//option == -1

      double a = M - 4.0;
      double b = (M - N - 2.0)/(N*(N - 1.0));
      double c = (M - N - 2.0)/(N - 1.0);

      this->Q(-1,a,b,c,tpm);

   }

}
