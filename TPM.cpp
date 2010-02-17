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

int TPM::counter = 0;

int **TPM::t2s;
int **TPM::s2t;

//constructor:
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

//copy constructor
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

//destructor
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

//access the numbers: sp indices
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

//friend function! output stream operator overloaded
ostream &operator<<(ostream &output,TPM &tpm_p){

   for(int i = 0;i < tpm_p.n;++i)
      for(int j = 0;j < tpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << tpm_p.t2s[i][0] << "\t" << tpm_p.t2s[i][1]

            << "\t" << tpm_p.t2s[j][0] << "\t" << tpm_p.t2s[j][1] << "\t" << tpm_p(i,j) << endl;

      }

   return output;

}

int TPM::gN(){

   return N;

}

int TPM::gM(){

   return M;

}

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

void TPM::Q(int option,TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

//een algemene Q-like afbeelding Q(A,B,C)(tpm_d)
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

void TPM::unit(){

   double ward = N*(N - 1.0)/(2.0*n);

   for(int i = 0;i < n;++i){

      (*this)(i,i) = ward;

      for(int j = i + 1;j < n;++j)
         (*this)(i,j) = (*this)(j,i) = 0.0;

   }

}

void TPM::proj_Tr(){

   double ward = (this->trace())/(double)n;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= ward;

}

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

//van de ph naar de tp ruimte:

//option == 1: doe dan de G down
//option == -1: doe dan de inverse G up
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

void TPM::min_unit(double scale){

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}

void TPM::min_qunit(double scale){

   double q = 1.0 + (M - 2*N)*(M - 1.0)/(N*(N - 1.0));

   scale *= q;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}
