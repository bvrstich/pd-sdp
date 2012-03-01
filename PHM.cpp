#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ifstream;
using std::endl;

#include "include.h"

int PHM::M;
int PHM::N;

vector< vector<int> > PHM::ph2s;
int **PHM::s2ph;

/**
 * initialize the static variables and allocate the static lists
 * @param M_i the number of sp orbitals
 * @param N_i the number of particles
 */
void PHM::init(int M_i,int N_i){

   M = M_i;
   N = N_i;

   //allocate
   s2ph = new int * [M];

   for(int a = 0;a < M;++a)
      s2ph[a] = new int [M];

   int ph = 0;

   vector<int> v(2);

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b){

         s2ph[a][b] = ph;

         v[0] = a;
         v[1] = b;

         ph2s.push_back(v);

         ++ph;

      }

}

/**
 * deallocate the lists
 */
void PHM::clear(){

   for(int a = 0;a < M;++a)
      delete [] s2ph[a];

   delete [] s2ph;

}

/**
 * standard constructor: constructs Matrix object of dimension M*M and
 */
PHM::PHM() : Matrix(ph2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*M and copies the content of phm_c into it,
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : Matrix(phm_c){ }

/**
 * destructor
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place PHM(i,j)
 */
double &PHM::operator()(int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);

}

/**
 * access the elements of the matrix in sp mode, 
 * @param a first sp index that forms the ph row index i together with b
 * @param b second sp index that forms the ph row index i together with a
 * @param c first sp index that forms the ph column index j together with d
 * @param d second sp index that forms the ph column index j together with c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int a,int b,int c,int d) const
{
   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(i,j);
}

ostream &operator<<(ostream &output,PHM &phm_p){

   for(int i = 0;i < phm_p.gn();++i)
      for(int j = 0;j < phm_p.gn();++j){

         output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

            << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return number of particles
 */
int PHM::gN() const
{
   return N;
}

/**
 * @return number of single particle oribals
 */
int PHM::gM() const
{
   return M;
}

/**
 * De G map, maps a TPM object on a PHM object.
 * @param option = 1 G_up map is used, = -1 G^{-1}_down map is used
 * @param tpm input TPM
 */
void PHM::G(int option,const TPM &tpm)
{
   SPM spm;

   if(option == 1)
      spm.bar(1.0/(N - 1.0),tpm);
   else
      spm.bar(1.0/(M - N + 1.0),tpm);

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < gn();++j){

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
 * The G2 map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G2(const TPM &tpm)
{
   SPM spm;
   spm.bar(1.0/(N - 1.0),tpm);

   double ward = tpm.trace()*2.0/(N*(N - 1.0));

   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < gn();++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = 0.0;

         if(a == b){

            if(c == d)
               (*this)(i,j) += ward;

            (*this)(i,j) -= spm(c,d);

         }

         if(c == d)
            (*this)(i,j) -= spm(a,b);

         if(a == c)
            (*this)(i,j) += spm(b,d);

         (*this)(i,j) -= tpm(a,d,c,b);

      }
   }
   
   this->symmetrize();

}

/**
 * Map a PPHM (pphm) object onto a PHM (*this) object by tracing one pair of indices (see primal_dual.pdf for more info)
 * @param pphm Input PPHM
 */
void PHM::bar(const PPHM &pphm)
{
   int a,b,c,d;

   for(int i = 0;i < gn();++i){

      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j = i;j < gn();++j){

         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < M;++l)
            (*this)(i,j) += pphm(l,a,b,l,c,d);

      }
   }

   this->symmetrize();

}

/**
 * fill the phm from a file with name filename, where the elements are indicated by their sp-indices
 * @param filename name of the inputfile
 */
void PHM::in_sp(const char *filename){

   ifstream input(filename);

   double value;

   int a,b,c,d;

   while(input >> a >> b >> c >> d >> value)
      (*this)(a,b,c,d) = value;

   this->symmetrize();

}

void PHM::bar(const T2PM &t2pm)
{
   int a,b,c,d;

   for(int i=0;i<gn();i++)
   {
      a = ph2s[i][0];
      b = ph2s[i][1];

      for(int j=i;j<gn();j++)
      {
         c = ph2s[j][0];
         d = ph2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l=0;l<M;l++)
            (*this)(i,j) += t2pm(l,a,b,l,c,d);

         (*this)(j,i) = (*this)(i,j);
      }
   }
}

/* vim: set ts=3 sw=3 expandtab :*/
