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

input *CartInt::readin;
vector< vector<int> > CartInt::s2inlxyz;
int ******CartInt::inlxyz2s;

int CartInt::l_max;
int CartInt::n_max;

vector< vector<int> > CartInt::t2s;
int **CartInt::s2t;

int CartInt::dim;

/** 
 * static function that reads in the input data and makes the matrix elements
 */
void CartInt::init(){

   readin = new input("start.stp");

   vector<int> v(6);

   int curtyp = 0;

   n_max = 0;
   l_max = 0;

   for(int i = 0;i < readin->gNcores();++i){

      v[0] = i;

      for(int j = 0;j < readin->gGaussInfo(i)->gNtypes();++j){

         if(readin->gGaussInfo(i)->gtype(j) == 'S')
            v[2] = 0;
         else if(readin->gGaussInfo(i)->gtype(j) == 'P')
            v[2] = 1;
         else if(readin->gGaussInfo(i)->gtype(j) == 'D')
            v[2] = 2;
         else if(readin->gGaussInfo(i)->gtype(j) == 'F')
            v[2] = 3;
         else if(readin->gGaussInfo(i)->gtype(j) == 'G')
            v[2] = 4;
         else if(readin->gGaussInfo(i)->gtype(j) == 'H')
            v[2] = 5;
         else
            cout << "BASISSET TOO LARGE" << endl;

         if(j == 0){

            v[1] = v[2] + 1;
            curtyp = v[2];

         }
         else{

            if(v[2] == curtyp)
               v[1]++;
            else{

               v[1] = v[2] + 1;
               curtyp = v[2];

            }

         }

         if(v[1] > n_max)
            n_max = v[1];

         if(v[2] > l_max)
            l_max = v[2];

         for(int x = v[2];x >= 0;x--)
            for(int y = v[2] - x;y >= 0;y--){

               v[3] = x;
               v[4] = y;
               v[5] = v[2] - x - y;

               s2inlxyz.push_back(v);

            }

      }

   }
   
   //allocate the list
   inlxyz2s =  new int ***** [readin->gNcores()];

   for(int i = 0;i < readin->gNcores();++i){

      inlxyz2s[i] = new int **** [n_max];

      for(int n = 0;n < n_max;++n){

         inlxyz2s[i][n] =  new int *** [l_max + 1];

         for(int l = 0;l <= l_max;++l){

            inlxyz2s[i][n][l] =  new int ** [l + 1];

            for(int x = 0;x <= l;++x){

               inlxyz2s[i][n][l][x] =  new int * [l + 1];

               for(int y = 0;y <= l;++y)
                  inlxyz2s[i][n][l][x][y] =  new int [l + 1];

            }

         }

      }

   }

   //fill list using other list
   for(unsigned int s = 0;s < s2inlxyz.size();++s){

      v = s2inlxyz[s];

      inlxyz2s[v[0]][v[1] - v[2] - 1][v[2]][v[3]][v[4]][v[5]] = s;

   }

   dim = s2inlxyz.size();

   s2t = new int * [dim];

   for(int i = 0;i < dim;++i)
      s2t[i] = new int [dim];

   vector<int> vst(2);

   int iter = 0;

   for(int i = 0;i < dim;++i)
      for(int j = 0;j < dim;++j){

         vst[0] = i;
         vst[1] = j;

         t2s.push_back(vst);

         s2t[i][j] = iter;

         ++iter;

      }

}

/** 
 * function that deallocates the static members
 */
void CartInt::clear(){

   for(int i = 0;i < readin->gNcores();++i){

      for(int n = 0;n < n_max;++n){

         for(int l = 0;l <= l_max;++l){

            for(int x = 0;x <= l;++x){

               for(int y = 0;y <= l;++y)
                  delete [] inlxyz2s[i][n][l][x][y];

               delete [] inlxyz2s[i][n][l][x];

            }

            delete [] inlxyz2s[i][n][l];

         }

         delete [] inlxyz2s[i][n];

      }

      delete [] inlxyz2s[i];

   }

   delete [] inlxyz2s;

   for(int i = 0;i < dim;++i)
      delete [] s2t[i];

   delete [] s2t;

   delete readin;

}

/** 
 * Standard constructor: allocates the matrices and fills them with the correct elements
 */
CartInt::CartInt(){ 

   S = new Matrix(dim);
   T = new Matrix(dim);
   U = new Matrix(dim);

   V = new Matrix(dim*dim);

   MxElem setup(*readin);
   setup.Init(*readin);

   for(int i = 0;i < dim;++i)
      for(int j = i;j < dim;++j){

         (*S)(i,j) = setup.gSoverlap(i,j);
         (*T)(i,j) = setup.gKEseparate(i,j);
         (*U)(i,j) = setup.gTelem(i,j) - setup.gKEseparate(i,j);

      }

   S->symmetrize();
   T->symmetrize();
   U->symmetrize();

   int a,b,c,d;

   for(int i = 0;i < dim*dim;++i){

      a = t2s[i][0];
      b = t2s[i][1];

      for(int j = i;j < dim*dim;++j){

         c = t2s[j][0];
         d = t2s[j][1];

         (*V)(i,j) = setup.gVelem(a,b,c,d);

      }
   }

   V->symmetrize();
   
}

/** 
 * copy constructor
 * @param ci_c CartInt object to be copied in the newly constructed object
 */
CartInt::CartInt(const CartInt &ci_c){ 

   S = new Matrix(ci_c.gS());
   T = new Matrix(ci_c.gT());
   U = new Matrix(ci_c.gU());
   
   V = new Matrix(ci_c.gV());

}

/**
 * standard destructor
 */
CartInt::~CartInt(){ 

   delete S;
   delete T;
   delete U;

   delete V;
   
}

/** 
 * @return the overlapmatrix, const version
 */
const Matrix &CartInt::gS() const { 

   return *S;

}

/** 
 * @return the overlapmatrix
 */
Matrix &CartInt::gS() { 

   return *S;

}

/** 
 * @return the kinetic energy matrix, const version
 */
const Matrix &CartInt::gT() const { 

   return *T; 
}

/** 
 * @return the kinetic energy matrix
 */
Matrix &CartInt::gT() { 

   return *T;

}

/** 
 * @return the nuclear attraction matrix, const version
 */
const Matrix &CartInt::gU() const { 

   return *U; 
}

/** 
 * @return the nuclear attraction matrix
 */
Matrix &CartInt::gU() { 

   return *U;

}

/** 
 * @return the electronic repulsion matrix
 */
const Matrix &CartInt::gV() const { 

   return *V; 
}

/** 
 * @return the electronic repulsion matrix
 */
Matrix &CartInt::gV() { 

   return *V;

}
