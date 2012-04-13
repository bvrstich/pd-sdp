#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;

#include "include.h"

vector< vector<int> > SphInt::s2inlm;
int ****SphInt::inlm2s;

vector< vector<int> > SphInt::t2s;
int **SphInt::s2t;

int SphInt::dim;
int SphInt::N;
int SphInt::N_Z;
int SphInt::n_max;
int SphInt::l_max;

/** 
 * static function that allocates the static lists and calculates the dimensions and such
 */
void SphInt::init(){

   N_Z = CartInt::gN_Z();
   n_max = CartInt::gn_max();
   l_max = CartInt::gl_max();
   N = CartInt::gN();

   //allocate
   inlm2s = new int *** [N_Z];

   for(int i = 0;i < N_Z;++i){

      inlm2s[i] = new int ** [n_max];

      for(int n = 0;n < n_max;++n){

         inlm2s[i][n] = new int * [l_max + 1];

         for(int l = 0;l <= l_max;++l)
            inlm2s[i][n][l] = new int [2*l + 1];

      }
   }

   //construct
   vector<int> v(4);

   for(int s = 0;s < CartInt::gdim();++s){

      v[0] = CartInt::gs2inlxyz(s,0);//i
      v[1] = CartInt::gs2inlxyz(s,1);//n
      v[2] = CartInt::gs2inlxyz(s,2);//l

      for(int m = -v[2];m <= v[2];++m){

         v[3] = m;//m

         s2inlm.push_back(v);

      }

      s += (v[2] + 2)*(v[2] + 1)/2 - 1;

   }

   dim = s2inlm.size();

   for(int s = 0;s < dim;++s){

      v = s2inlm[s];

      inlm2s[v[0]][v[1] - v[2] - 1][v[2]][v[3] + v[2]] = s;

   }

}

/** 
 * function that deallocates the static members
 */
void SphInt::clear(){

   for(int i = 0;i < N_Z;++i){

      for(int n = 0;n < n_max;++n){

         for(int l = 0;l <= l_max;++l)
            delete [] inlm2s[i][n][l];

         delete [] inlm2s[i][n];

      }

      delete [] inlm2s[i];

   }

   delete [] inlm2s;
    
}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements by transforming a CartInt object
 * @param ci input CartInt object
 */
SphInt::SphInt(const CartInt &ci){ 

   S = new Matrix(dim);
   T = new Matrix(dim);
   U = new Matrix(dim);

   V = new Matrix(dim*dim);

   int i,n_i,l_i,m_i;
   int j,n_j,l_j,m_j;

   //start with overlap
   for(int s_i = 0;s_i < dim;++s_i){

      i = s2inlm[s_i][0];
      n_i = s2inlm[s_i][1];
      l_i = s2inlm[s_i][2];
      m_i = s2inlm[s_i][3];

      for(int s_j = s_i;s_j < dim;++s_j){

         j = s2inlm[s_j][0];
         n_j = s2inlm[s_j][1];
         l_j = s2inlm[s_j][2];
         m_j = s2inlm[s_j][3];

         (*S)(s_i,s_j) = 0.0;

         if(i == j){//basisfunctions on the same core

            if(l_i == l_j && m_i == m_j){

               if(l_i == 0){

                  (*S)(s_i,s_j) = ci.gS(i,n_i,l_i,0,0,0,j,n_j,l_j,0,0,0);

               }
               else if(l_i == 1){

                  if(m_i == 0)
                     (*S)(s_i,s_j) = ci.gS(i,n_i,l_i,0,0,1,j,n_j,l_j,0,0,1);
                  else
                     (*S)(s_i,s_j) = 0.5 * ( ci.gS(i,n_i,l_i,1,0,0,j,n_j,l_j,1,0,0) + ci.gS(i,n_i,l_i,0,1,0,j,n_j,l_j,0,1,0) );

               }
               else if(l_i == 2){

                  if(m_i == 0){

                     (*S)(s_i,s_j) = ci.gS(i,n_i,l_i,0,0,2,j,n_j,l_j,0,0,2)
                     
                        + 0.25 * ( ci.gS(i,n_i,l_i,0,2,0,j,n_j,l_j,0,2,0) +  ci.gS(i,n_i,l_i,2,0,0,j,n_j,l_j,2,0,0) 
                        
                              + 2.0 * ci.gS(i,n_i,l_i,0,2,0,j,n_j,l_j,2,0,0) )

                        - ci.gS(i,n_i,l_i,0,0,2,j,n_j,l_j,2,0,0) - ci.gS(i,n_i,l_i,0,0,2,j,n_j,l_j,0,2,0);

                  }
                  else if(m_i == 1 || m_i == -1){

                     (*S)(s_i,s_j) = 0.5 * ( ci.gS(i,n_i,l_i,1,0,1,j,n_j,l_j,1,0,1) + ci.gS(i,n_i,l_i,0,1,1,j,n_j,l_j,0,1,1) );

                  }
                  else if(m_i == 2 || m_i == -2){

                     (*S)(s_i,s_j) = 3.0/8.0 * ( ci.gS(i,n_i,l_i,2,0,0,j,n_j,l_j,2,0,0) + ci.gS(i,n_i,l_i,0,2,0,j,n_j,l_j,0,2,0)
                     
                           - 2.0 * ci.gS(i,n_i,l_i,2,0,0,j,n_j,l_j,0,2,0)) + 0.5 * ci.gS(i,n_i,l_i,1,1,0,j,n_j,l_j,1,1,0);

                  }

               }
               else if(l_i == 3){

                  cout << "If I get here, add more" << endl;

               }
               else if(l_i == 4){

                  cout << "If I get here, add more" << endl;

               }
               else if(l_i == 5){

                  cout << "If I get here, add more" << endl;

               }
               else
                  cout << "Basisset too large for me" << endl;

            }

         }
         else{//basisfunctions on different cores

         }

      }
   }

}

/** 
 * copy constructor
 * @param ci_c SphInt object to be copied in the newly constructed object
 */
SphInt::SphInt(const SphInt &ci_c){ 

   S = new Matrix(ci_c.gS());
   T = new Matrix(ci_c.gT());
   U = new Matrix(ci_c.gU());

   V = new Matrix(ci_c.gV());

}

/**
 * standard destructor
 */
SphInt::~SphInt(){ 

   delete S;
   delete T;
   delete U;

   delete V;

}

/** 
 * @return the overlapmatrix, const version
 */
const Matrix &SphInt::gS() const { 

   return *S;

}

/** 
 * @return the overlapmatrix
 */
Matrix &SphInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix, const version
 */
const Matrix &SphInt::gT() const { 

   return *T; 
}

/** 
 * @return the kinetic energy matrix
 */
Matrix &SphInt::gT() { 

   return *T;

}

/** 
 * @return the nuclear attraction matrix, const version
 */
const Matrix &SphInt::gU() const { 

   return *U; 
}

/** 
 * @return the nuclear attraction matrix
 */
Matrix &SphInt::gU() { 

   return *U;

}

/** 
 * @return the electronic repulsion matrix
 */
const Matrix &SphInt::gV() const { 

   return *V; 
}

/** 
 * @return the electronic repulsion matrix
 */
Matrix &SphInt::gV() { 

   return *V;

}

/**
 * @return the dimension of spatial sp space
 */
int SphInt::gdim() {

   return dim;

}

/**
 * static function
 * @return nr of electrons
 */
int SphInt::gN(){

   return N;

}
