#include <iostream>

#include "Math.h"

/**
 * just a simple assisting function for calculating combinations.
 * @return the number of combinations of N out of M.
 */
int Math::comb(int M,int N){

   int max;
   int min;

   if(N > M/2){

      max = N;
      min = M - N;

   }
   else{

      max = M - N;
      min = N;

   }

   int tmp = M;

   for(int i = M - 1;i > max;--i)
      tmp *= i;

   for(int i = min;i > 1;--i)
      tmp /= i;

   return tmp;

}

/**
 * Function that calculates the adjoint sp-orbital, so turns spin up into spin down on same spatial orbital.
 * @param index the original index
 * @return the adjoint index
 */
int Math::adjoint(int index){

   if(index % 2 == 0)
      return index + 1;
   else
      return index - 1;

}
