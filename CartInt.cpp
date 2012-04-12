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

   delete readin;

}

/** 
 * Standard constructor
 */
CartInt::CartInt(){ }

/**
 * standard destructor
 */
CartInt::~CartInt(){ }
