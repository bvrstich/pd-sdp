#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int T2PM::counter = 0;

int **T2PM::pph2s;
int ***T2PM::s2pph;

/**
 * standard constructor: constructs Matrix object of dimension M+M^2*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and T2' basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
T2PM::T2PM(int M,int N) : Matrix(M+M*M*(M-1)/2)
{
   this->N = N;
   this->M = M;
   this->n_pph = M*M*(M-1)/2;
   this->n = M + n_pph;

   if(counter == 0)
   {
      //allocatie van s2pph
      s2pph = new int ** [M];

      for(int i = 0;i < M;i++)
      {
         s2pph[i] = new int * [M];

         for(int j = 0;j < M;j++)
            s2pph[i][j] = new int [M];
      }

      //allocatie van pph2s
      pph2s = new int * [n];

      for(int i = 0;i < n;i++)
         pph2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;a++)
         for(int b = a + 1;b < M;b++)
            for(int c = 0;c < M;c++)
            {
               s2pph[a][b][c] = teller;

               pph2s[teller][0] = a;
               pph2s[teller][1] = b;
               pph2s[teller++][2] = c;
            }
   }
   counter++;
}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)*(M - 2)/6 and copies the content of T2PM_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param T2PM_c input T2PM to be copied
 */
T2PM::T2PM(const T2PM &T2PM_c) : Matrix(T2PM_c)
{
   this->N = T2PM_c.N;
   this->M = T2PM_c.M;
   this->n = T2PM_c.n;
   this->n_pph = T2PM_c.n_pph;

   if(counter == 0)
   {
      //allocatie van s2pph
      s2pph = new int ** [M];

      for(int i = 0;i < M;i++)
      {
         s2pph[i] = new int * [M];

         for(int j = 0;j < M;j++)
            s2pph[i][j] = new int [M];
      }

      //allocatie van pph2s
      pph2s = new int * [n];

      for(int i = 0;i < n;++i)
         pph2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = a + 1;b < M;++b)
            for(int c = 0;c < M;++c)
            {
               s2pph[a][b][c] = teller;

               pph2s[teller][0] = a;
               pph2s[teller][1] = b;
               pph2s[teller++][2] = c;
            }
   }
   counter++;
}

/**
 * destructor: if counter == 1 the memory for the static lists pph2s en s2pph twill be deleted.
 */
T2PM::~T2PM()
{
   if(counter == 1)
   {
      for(int i = 0;i < M;++i)
      {
         for(int j = 0;j < M;++j)
            delete [] s2pph[i][j];

         delete [] s2pph[i];
      }

      delete [] s2pph;

      for(int i = 0;i < n;++i)
         delete [] pph2s[i];

      delete [] pph2s;
   }

   --counter;
}

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
/*    output << std::setprecision(2) << std::fixed;
 * 
 *    for(int i = 0;i < T2PM_p.n;i++)
 *    {
 *       for(int j = 0;j < T2PM_p.n;j++)
 *          output << std::setfill('0') << std::setw(6) << T2PM_p(i,j) << " ";
 * 
 *       output << endl;
 *    }
 * 
 *    output << endl;
 */
   output << std::setprecision(10) << std::scientific;

   // T2 part
   for(int i = 0;i < T2PM_p.n_pph;i++)
      for(int j = i;j <  T2PM_p.n_pph;j++)
         output << i << "\t" << j << "\t|\t" << T2PM_p.pph2s[i][0] << "\t" << T2PM_p.pph2s[i][1] << "\t"
            << T2PM_p.pph2s[i][2] << "\t" << T2PM_p.pph2s[j][0] << "\t" << T2PM_p.pph2s[j][1] << "\t"
            << T2PM_p.pph2s[j][2] << "\t -> " << T2PM_p(i,j) << endl;

   // off diagonal part
   for(int i = 0;i < T2PM_p.n_pph;i++)
      for(int j = 0;j <  T2PM_p.M;j++)
         output << i << "\t" << j << "\t|\t" << T2PM_p.pph2s[i][0] << "\t" << T2PM_p.pph2s[i][1] << "\t"
            << T2PM_p.pph2s[i][2] << "\t" << j << "\t -> " << T2PM_p(i,j+T2PM_p.n_pph)
            << endl;

   // rho part
   for(int i = 0;i < T2PM_p.M;i++)
      for(int j = i;j <  T2PM_p.M;j++)
         output << i+T2PM_p.n_pph << "\t" << j+T2PM_p.n_pph << "\t|\t" << i << "\t" << j << "\t -> " << T2PM_p(i+T2PM_p.n_pph,j+T2PM_p.n_pph)
            << endl;

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
 * @return dimension of T2' space and of Matrix
 */
int T2PM::gn() const
{
   return n;
}

/**
 * Deduct scale times the T2P of the unit matrix from (*this).
 * @param scale The number by which to scale the unitmatrix.
 */
void T2PM::min_tunit(double scale)
{
   int i,j;

   // T2, the off-diagonal part
   for(int a = 0;a < M;++a)
   {
      //first a > b
      for(int b = 0;b < a;++b)
         for(int c = a;c < M;++c)
         {
            //c always >= a
            i = s2pph[b][a][a];
            j = s2pph[b][c][c];

            (*this)(i,j) -= scale;
         }

      //then a < b
      for(int b = a + 1;b < M;++b)
      {
         //first c < b
         for(int c = a;c < b;++c)
         {
            i = s2pph[a][b][a];
            j = s2pph[c][b][c];

            (*this)(i,j) -= scale;
         }

         //then c > b
         for(int c = b + 1;c < M;++c)
         {
            i = s2pph[a][b][a];
            j = s2pph[b][c][c];

            (*this)(i,j) += scale;
         }
      }
   }

   // T2, the diagonal part
   double t2p = (M - N)/(N - 1.0) * scale;

   for(i=0;i<n_pph;i++)
      (*this)(i,i) -= t2p;

   // the rho part
   t2p = (M-1.0)/(N-1.0) * scale;

   for(i=n_pph;i<n;i++)
      (*this)(i,i) -= t2p;

   // the omega part
   for(i=0;i<n_pph;i++)
      for(j=0;j<M;j++)
      {
         int a,b,c,d;
         a = pph2s[i][0];
         b = pph2s[i][1];
         c = pph2s[i][2];
         d = j;

         if( a==d && b==c )
            (*this)(i,n_pph+j) -= scale;

         if( a==c && b==d )
            (*this)(i,n_pph+j) += scale;
      }

   this->symmetrize();
}

/** 
 * @return the skew trace, for T2PM matrices defined as sum_abc T2PM(a,b,a,c,b,c)
 */
double T2PM::skew_trace() const
{
   double ward = 0.0;

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b)
         for(int c = 0;c < M;++c)
            ward += (*this)(a,b,a,c,b,c);

   return ward;
}

/**
 * @return the trace over the T2 part
 */
double T2PM::T2_trace() const
{
   double brecht = 0;

   for(int i=0;i<n_pph;i++)
      brecht += (*this)(i,i);

   return brecht;
}

/**
 * @return the trace over the rho part
 */
double T2PM::rho_trace() const
{
   double brecht = 0;

   for(int i=n_pph;i<n;i++)
      brecht += (*this)(i,i);

   return brecht;
}

/**
 * @return the skew trace over the omega part, defined as sum_abc T2PM(a,b,b,a)-T2PM(a,b,a,b)
 */
double T2PM::omega_trace() const
{
   double brecht = 0;

   for(int a=0;a<M;a++)
      for(int b=a+1;b<M;b++)
         brecht += (*this)(a,b,b,a) - (*this)(a,b,a,b);

   return brecht;
}

/**
 * The T2-map: maps a TPM object (tpm) on a T2PM object (*this)
 * @param tpm input TPM
 */
void T2PM::T(const TPM &tpm)
{
   int a,b,c,d,e,z;

   SPM spm(1.0/(N-1.0),tpm);

   for(int i = 0;i < n_pph;i++)
   {
      a = pph2s[i][0];
      b = pph2s[i][1];
      c = pph2s[i][2];

      for(int j = i;j < n_pph;j++)
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

         (*this)(n_pph+i,n_pph+j) = spm(i,j);

         (*this)(n_pph+j,n_pph+i) = (*this)(n_pph+i,n_pph+j);

      }


   for(int i = 0;i < n_pph;i++)
   {
      a = pph2s[i][0];
      b = pph2s[i][1];
      c = pph2s[i][2];

      for(int j = 0;j < M;j++)
      {
         (*this)(i,n_pph+j) = tpm(a,b,j,c);

         (*this)(n_pph+j,i) = (*this)(i,j+n_pph);
      }
   }
}

/**
 * Cast operator: this implements a cast from a T2PM object to a PPHM object by copying the A_1 matrix of T2PM to PPHM
 */
T2PM::operator PPHM() const
{
   PPHM A(M,N);

   for(int i = 0;i < n_pph;i++)
      for(int j = i;j < n_pph;j++)
         A(i,j) = A(j,i) = (*this)(i,j);

   return A;
}

/**
 * Cast operator: this implements a cast from a T2PM object to a SPM object by copying the A_2 matrix of T2PM to SPM
 */
T2PM::operator SPM() const
{
   SPM A(M,N);

   for(int i = 0;i < M;i++)
      for(int j = i;j < M;j++)
         A(i,j) = A(j,i) = (*this)(n_pph+i,n_pph+j);

   return A;
}

/* vim: set ts=3 sw=3 expandtab :*/
