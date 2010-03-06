#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int PHM::counter = 0;

int **PHM::ph2s;
int **PHM::s2ph;

/**
 * standard constructor, maakt Matrix aan met dimensie M*M
 * als counter == 0 dan alloceerd de constructor het geheugen voor lijsten ph2s en s2ph en initialiseerd deze:
 * @param M aantal sp orbitals
 * @param N aantal deeltjes
 */
PHM::PHM(int M,int N) : Matrix(M*M) {
   
   this->N = N;
   this->M = M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van s2ph
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van ph2s
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * copy constructor, maakt Matrix aan met dimensie M*M en kopieerd er phm_c in
 * als counter == 0 dan alloceerd de constructor het geheugen voor lijsten ph2s en s2ph en initialiseerd deze:
 * @param phm_c matrix die gekopieerd moet worden
 */
PHM::PHM(PHM &phm_c) : Matrix(phm_c){

   this->N = phm_c.N;
   this->M = phm_c.M;
   this->n = M*M;

   if(counter == 0){

      //allocatie van sp2tp
      s2ph = new int * [M];
      s2ph[0] = new int [M*M];

      for(int i = 1;i < M;++i)
         s2ph[i] = s2ph[i - 1] + M;

      //allocatie van tp2sp
      ph2s = new int * [n];

      for(int i = 0;i < n;++i)
         ph2s[i] = new int [2];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b){

            s2ph[a][b] = teller;

            ph2s[teller][0] = a;
            ph2s[teller][1] = b;

            ++teller;

         }

   }

   ++counter;

}

/**
 * Destructor: als counter == 1 dan dealloceerd de destructor het geheugen voor lijsten ph2s en s2ph.
 */
PHM::~PHM(){

   if(counter == 1){

      delete [] s2ph[0];
      delete [] s2ph;

      for(int i = 0;i < n;++i)
         delete [] ph2s[i];

      delete [] ph2s;

   }

   --counter;

}

/**
 * toegang tot de getallen in de PHM door gebruik te maken van de sp indices,
 * @param a eerste sp index, vormt samen met b de eerste ph index
 * @param b tweede sp index, vormt samen met a de eerste ph index
 * @param c derde sp index, vormt samen met d de tweede ph index
 * @param d vierde sp index, vormt samen met c de tweede ph index
 */
double &PHM::operator()(int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,PHM &phm_p){

   for(int i = 0;i < phm_p.n;++i)
      for(int j = 0;j < phm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

            << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return aantal deeltjes
 */
int PHM::gN(){

   return N;

}

/**
 * @return aantal sp orbitals
 */
int PHM::gM(){

   return M;

}

/**
 * @return de dimensie van de particle hole ruimte en tevens de dimensie van de matrix
 */
int PHM::gn(){

   return n;

}

/**
 * De G afbeelding die een TPM object afbeeld op een PHM object.
 * @param option = 1 dan wordt G_up uitgevoerd, = -1 dan wordt G^{-1}_down uitgevoerd
 * @param tpm input TPM die afgebeeld wordt op this
 */
void PHM::G(int option,TPM &tpm){

   SPM spm(M,N);

   if(option == 1)
      spm.constr(1.0/(N - 1.0),tpm);
   else
      spm.constr(1.0/(M - N + 1.0),tpm);

   int a,b,c,d;

   for(int i = 0;i < n;++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < n;++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = -tpm(a,d,c,b);

         if(b == d)
            (*this)(i,j) += spm(a,c);

      }
   }
   
   //nog schalen met 4 voor G^{-1}, door de G down die eigenlijk een factor 4 te groot is
   if(option == -1)
      this->dscal(0.25);

   this->symmetrize();

}

/**
 * Bereken de skew trace gedefinieerd als: \n
 * sum_{ab} PHM(a,a,b,b)
 * @return de skew trace
 */
double PHM::skew_trace(){

   double ward = 0.0;

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         ward += (*this)(a,a,b,b);

   return ward;

}

/**
 * Trek van this - de G-afbeelding van de eenheidsmatrix  * een constante - af:\n
 * this -= scale* G(1) \n
 * zie ook de nota's primal_dual.pdf voor meer info.
 * @param scale de constante waarmee de eenheidsmatrix vermenigvuldigt wordt
 */
void PHM::min_gunit(double scale){

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         (*this)(a,a,b,b) -= scale;

   double g = (M - N)/(N - 1.0);

   scale *= g;

   for(int i = 0;i < n;++i)
      (*this)(i,i) -= scale;

}
