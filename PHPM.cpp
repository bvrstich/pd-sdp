#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int PHPM::counter = 0;

int **PHPM::php2s;
int ***PHPM::s2php;

/**
 * standard constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and php basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PHPM::PHPM(int M,int N) : Matrix(M*M*M) {

   this->N = N;
   this->M = M;
   this->n = M*M*M;

   if(counter == 0){

      //allocatie van s2php
      s2php = new int ** [M];

      for(int i = 0;i < M;++i){

         s2php[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2php[i][j] = new int [M];

      }

      //allocatie van php2s
      php2s = new int * [n];

      for(int i = 0;i < n;++i)
         php2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               s2php[a][b][c] = teller;

               php2s[teller][0] = a;
               php2s[teller][1] = b;
               php2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*M*(M - 1)/2 and copies the content of phpm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and php basis.
 * @param phpm_c input PHPM to be copied
 */
PHPM::PHPM(const PHPM &phpm_c) : Matrix(phpm_c){

   this->N = phpm_c.N;
   this->M = phpm_c.M;
   this->n = phpm_c.n;

   if(counter == 0){

      //allocatie van s2php
      s2php = new int ** [M];

      for(int i = 0;i < M;++i){

         s2php[i] = new int * [M];

         for(int j = 0;j < M;++j)
            s2php[i][j] = new int [M];

      }

      //allocatie van php2s
      php2s = new int * [n];

      for(int i = 0;i < n;++i)
         php2s[i] = new int [3];

      //initialisatie van de twee arrays
      int teller = 0;

      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c){

               s2php[a][b][c] = teller;

               php2s[teller][0] = a;
               php2s[teller][1] = b;
               php2s[teller][2] = c;

               ++teller;

            }

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists php2s en s2php will be deleted.
 */
PHPM::~PHPM(){

   if(counter == 1){

      for(int i = 0;i < M;++i){

         for(int j = 0;j < M;++j)
            delete [] s2php[i][j];

         delete [] s2php[i];

      }

      delete [] s2php;

      for(int i = 0;i < n;++i)
         delete [] php2s[i];

      delete [] php2s;

   }

   --counter;

}

/**
 * access the elements of the matrix in sp mode
 * PHPM(a,b,c,d,e,f) = -PHPM(b,a,c,d,e,f)  but PHPM(a,b,c,d,e,f) != - PHPM(a,c,b,d,e,f) \n\n
 * PHPM(a,a,c,d,e,f) = 0 but PHPM(a,b,b,d,e,f) != 0 \n\n
 * @param a first sp index that forms the php row index i together with b and c
 * @param b second sp index that forms the php row index i together with a and c
 * @param c third sp index that forms the php row index i together with a and b
 * @param d first sp index that forms the php column index j together with e and z
 * @param e second sp index that forms the php column index j together with d and z
 * @param z third sp index that forms the php column index j together with d and e
 * @return the number on place PHPM(i,j) with the right phase.
 */
double PHPM::operator()(int a,int b,int c,int d,int e,int z) const{

   return (*this)(s2php[a][b][c],s2php[d][e][z]);

}

ostream &operator<<(ostream &output,PHPM &phpm_p){

   for(int i = 0;i < phpm_p.n;++i)
      for(int j = 0;j < phpm_p.n;++j){

         output << i << "\t" << j << "\t|\t" << phpm_p.php2s[i][0] << "\t" << phpm_p.php2s[i][1] << "\t" << phpm_p.php2s[i][2]

            << "\t" << phpm_p.php2s[j][0] << "\t" << phpm_p.php2s[j][1] << "\t" << phpm_p.php2s[j][2] << "\t" << phpm_p(i,j) << endl;

      }

   return output;

}

/**
 * @return nr of particles
 */
int PHPM::gN() const
{
   return N;
}

/**
 * @return dimension of sp space
 */
int PHPM::gM() const
{
   return M;
}

/**
 * @return dimension of php space and of Matrix
 */
int PHPM::gn() const
{
   return n;
}

/**
 * The T3-map: maps a TPM object (tpm) on a PHPM object (*this). see higher_order.pdf for more information
 * @param tpm input TPM
 */
void PHPM::T(const TPM &tpm){

   double ward = 2.0 * tpm.trace()/ (N * (N - 1.0));

   //construct the spm
   SPM spm(1.0/(N - 1.0),tpm);

   int a,b,c,d,e,z;

   for(int i = 0;i < n;++i){

      a = php2s[i][0];
      b = php2s[i][1];
      c = php2s[i][2];

      for(int j = i;j < n;++j){

         d = php2s[j][0];
         e = php2s[j][1];
         z = php2s[j][2];

         //initialize
         (*this)(i,j) = 0.0;

         if(a == d){

            (*this)(i,j) -= tpm(e,c,b,z);

            if(e == z)
               (*this)(i,j) -= spm(b,c);

            if(c == z)
               (*this)(i,j) += spm(b,e);

            if(b == c){

               (*this)(i,j) -= spm(e,z);

               if(e == z)
                  (*this)(i,j) += ward;

            }

         }

         if(a == z){

            (*this)(i,j) -= tpm(e,c,d,b);

            if(b == c)
               (*this)(i,j) += spm(e,d);

            if(c == d)
               (*this)(i,j) -= spm(b,e);

         }

         if(c == d){

            (*this)(i,j) -= tpm(a,e,b,z);

            if(e == z)
               (*this)(i,j) += spm(a,b);

         }

         if(b == e)
            (*this)(i,j) += tpm(a,c,d,z);

         if(c == z)
            (*this)(i,j) -= tpm(a,e,d,b);

      }
   }

   //and symmetrize
   this->symmetrize();

}

/* vim: set ts=3 sw=3 expandtab :*/
