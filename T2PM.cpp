#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

vector< vector<int> > T2PM::pph2s;
int ***T2PM::s2pph;

int T2PM::M;
int T2PM::N;

/**
 * initialize the static variables and allocate the static lists
 * @param M_i the number of sp orbitals
 * @param N_i the number of particles
 */
void T2PM::init(int M_i,int N_i){

   M = M_i;
   N = N_i;

   //allocate
   s2pph = new int ** [M];

   for(int a = 0;a < M;++a){

      s2pph[a] = new int * [M];

      for(int b = 0;b < M;++b)
         s2pph[a][b] = new int [M];

   }

   int pph = 0;

   vector<int> v(3);

   for(int a = 0;a < M;++a)
      for(int b = a + 1;b < M;++b)
         for(int c = 0;c < M;++c){

            s2pph[a][b][c] = pph;

            v[0] = a;
            v[1] = b;
            v[2] = c;

            pph2s.push_back(v);

            ++pph;

         }

}

/**
 * deallocate the static lists
 */
void T2PM::clear(){

   for(int a = 0;a < M;++a){

      for(int b = 0;b < M;++b)
         delete [] s2pph[a][b];

      delete [] s2pph[a];

   }

   delete [] s2pph;

}

/**
 * standard constructor: constructs Matrix object of dimension M+M^2*(M - 1)/2 and
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
T2PM::T2PM() : Matrix(pph2s.size() + M) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)*(M - 2)/6 and copies the content of T2PM_c into it,
 * @param T2PM_c input T2PM to be copied
 */
T2PM::T2PM(const T2PM &T2PM_c) : Matrix(T2PM_c) { }

/**
 * destructor
 */
T2PM::~T2PM() { }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * T2PM(a,b,c,d,e,f) = -T2PM(b,a,c,d,e,f) = ...\n\n
 * T2PM(a,a,c,d,e,f) = 0\n\n
 * @param a first sp index that forms the pph row index i together with b and c
 * @param b second sp index that forms the pph row index i together with a and c
 * @param c third sp index that forms the pph row index i together with a and b
 * @param d first sp index that forms the pph column index j together with e and z
 * @param e second sp index that forms the pph column index j together with d and z
 * @param z third sp index that forms the pph column index j together with d and e
 * @return the number on place T2PM(i,j) with the right phase.
 */
double T2PM::operator()(int a,int b,int c,int d,int e,int z) const
{
   //eerst kijken of er geen indices gelijk zijn:
   if(a == b || d == e)
      return 0;

   //dan kijken wel pph index met welke fase moet genomen worden:
   int i,j;

   int phase = 1;

   if(a < b)
      i = s2pph[a][b][c];
   else
   {
      i = s2pph[b][a][c];
      phase *= -1;
   }

   if(d < e)
      j = s2pph[d][e][z];
   else
   {
      j = s2pph[e][d][z];
      phase *= -1;
   }

   return phase*(*this)(i,j);
}

double T2PM::operator()(int a,int b,int c,int d) const
{
   //eerst kijken of er geen indices gelijk zijn:
   if(a == b)
      return 0;

   //dan kijken wel pph index met welke fase moet genomen worden:
   int i;

   int phase = 1;

   if(a < b)
      i = s2pph[a][b][c];
   else
   {
      i = s2pph[b][a][c];
      phase *= -1;
   }

   return phase*(*this)(i,d+M*M*(M-1)/2);
}

ostream &operator<<(ostream &output,const T2PM &T2PM_p)
{

   output << std::setprecision(10) << std::scientific;

   // T2 part
   for(unsigned int i = 0;i < T2PM_p.pph2s.size();i++)
      for(unsigned int j = i;j <  T2PM_p.pph2s.size();j++)
         output << i << "\t" << j << "\t|\t" << T2PM_p.pph2s[i][0] << "\t" << T2PM_p.pph2s[i][1] << "\t"

            << T2PM_p.pph2s[i][2] << "\t" << T2PM_p.pph2s[j][0] << "\t" << T2PM_p.pph2s[j][1] << "\t"

            << T2PM_p.pph2s[j][2] << "\t -> " << T2PM_p(i,j) << endl;

   // off diagonal part
   for(unsigned int i = 0;i < T2PM_p.pph2s.size();i++)
      for(int j = 0;j < T2PM_p.gM();j++)
         output << i << "\t" << j << "\t|\t" << T2PM_p.pph2s[i][0] << "\t" << T2PM_p.pph2s[i][1] << "\t"

            << T2PM_p.pph2s[i][2] << "\t" << j << "\t -> " << T2PM_p(i,j+T2PM_p.pph2s.size()) << endl;

   // rho part
   for(int i = 0;i < T2PM_p.gM();i++)
      for(int j = i;j <  T2PM_p.gM();j++)
         output << i+T2PM_p.pph2s.size() << "\t" << j+T2PM_p.pph2s.size() << "\t|\t" << i << "\t" << j << "\t -> " 
         
         << T2PM_p(i+T2PM_p.pph2s.size(),j+T2PM_p.pph2s.size()) << endl;

   output.unsetf(std::ios_base::floatfield);

   return output;
}

/**
 * @return nr of particles
 */
int T2PM::gN() const
{
   return N;
}

/**
 * @return dimension of sp space
 */
int T2PM::gM() const
{
   return M;
}

/**
 * The T2-map: maps a TPM object (tpm) on a T2PM object (*this)
 * @param tpm input TPM
 */
void T2PM::T(const TPM &tpm)
{
   int a,b,c,d,e,z;

   SPM spm;
   spm.bar(1.0/(N-1.0),tpm);

   for(unsigned int i = 0;i < pph2s.size();i++)
   {
      a = pph2s[i][0];
      b = pph2s[i][1];
      c = pph2s[i][2];

      for(unsigned int j = i;j < pph2s.size();j++)
      {
         d = pph2s[j][0];
         e = pph2s[j][1];
         z = pph2s[j][2];

         (*this)(i,j) = 0;

         // sp terms
         if( a == d && b == e )
            (*this)(i,j) += spm(c,z);

         if( a == e && b == d )
            (*this)(i,j) -= spm(c,z);

         // tp terms
         if( a == d )
            (*this)(i,j) -= tpm(z,b,c,e);

         if( b == d )
            (*this)(i,j) += tpm(z,a,c,e);

         if( c == z )
            (*this)(i,j) += tpm(a,b,d,e);

         if( b == e )
            (*this)(i,j) -= tpm(z,a,c,d);

         if( a == e )
            (*this)(i,j) += tpm(z,b,c,d);


         (*this)(j,i) = (*this)(i,j);
      }
   }

   for(int i = 0;i < M;i++)
      for(int j = i;j < M;j++)
      {

         (*this)(pph2s.size()+i,pph2s.size()+j) = spm(i,j);

         (*this)(pph2s.size()+j,pph2s.size()+i) = (*this)(pph2s.size()+i,pph2s.size()+j);

      }


   for(unsigned int i = 0;i < pph2s.size();i++)
   {
      a = pph2s[i][0];
      b = pph2s[i][1];
      c = pph2s[i][2];

      for(int j = 0;j < M;j++)
      {
         (*this)(i,pph2s.size()+j) = tpm(a,b,j,c);

         (*this)(pph2s.size()+j,i) = (*this)(i,j+pph2s.size());
      }
   }
}

/**
 * Cast operator: this implements a cast from a T2PM object to a PPHM object by copying the A_1 matrix of T2PM to PPHM
 */
T2PM::operator PPHM() const
{
   PPHM A;

   for(unsigned int i = 0;i < pph2s.size();i++)
      for(unsigned int j = i;j < pph2s.size();j++)
         A(i,j) = A(j,i) = (*this)(i,j);

   return A;
}

/**
 * Cast operator: this implements a cast from a T2PM object to a SPM object by copying the A_2 matrix of T2PM to SPM
 */
T2PM::operator SPM() const
{
   SPM A;

   for(int i = 0;i < M;i++)
      for(int j = i;j < M;j++)
         A(i,j) = A(j,i) = (*this)(pph2s.size()+i,pph2s.size()+j);

   return A;
}

/* vim: set ts=3 sw=3 expandtab :*/
