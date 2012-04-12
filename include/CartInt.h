#ifndef CARTINT_H
#define CARTINT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "input.h"

/**
 * @author Brecht Verstichel
 * @date 12-04-2012\n\n
 * This class CartInt is a class written for the storage of matrixelements in a cartesian Gaussian basis
 * The goal is to translate ThING input to a framework where I can transform it to a spherical basis.
 */
class CartInt {

   /**
    * Output stream operator overloaded
    * @param cartint_p the CartInt you want to print
    */
   friend ostream &operator<<(ostream &output,CartInt &cartint_p);

   public:
      
      //constructor
      CartInt();

      //copy constructor
      CartInt(const CartInt &);

      //destructor
      virtual ~CartInt();

      static void init();

      static void clear();

   private:

      //!static objects needed to construct and destruct all the lists
      static int l_max,n_max;

      //!input object contains all info about system
      static input *readin;

      //!list to switch between matrix index and physical quantum numbers
      static vector< vector<int> > s2inlxyz;

      //!list to switch between matrix index and physical quantum numbers
      static int ******inlxyz2s;

      //!dimension of the basisset
      int dim;
      
      //!overlapmatrix
      Matrix *S;

      //!kinetic energy matrix
      Matrix *T;

      //!nuclear attraction energy matrix
      Matrix *U;

      //!electronic repulsion matrix
      Matrix *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
