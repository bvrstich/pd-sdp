#ifndef EIG_PQ_H
#define EIG_PQ_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "../SUP/SUP_PQ.h"

/**
 * @author Brecht Verstichel
 * @date 04-03-2010\n\n
 * This class, EIG_PQ is a "block"-vector over the carrierspace's of the P and Q conditions. It contains room
 * to store the eigenvalues of a P and a Q block, and special member function that work with these eigenvalues.
 * This class should only be used when a SUP_PQ matrix has been diagonalized, some functions could give strange results when the EIG_PQ object is filled
 * with random numbers.\n\n
 * This function is the mother class of all EIG_PQ* classes.
 */
class EIG_PQ{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG_PQ you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG_PQ &eig_p);

   public:
      
      //constructor
      EIG_PQ(int M,int N);

      //copy constructor
      EIG_PQ(EIG_PQ &);

      //constructor met initialisatie op 
      EIG_PQ(SUP_PQ &);

      //destructor
      ~EIG_PQ();

      int gN();

      int gM();

      int gn_tp();

      double centerpot(double,EIG_PQ &,double,double);

      double operator()(int,int);

      //overload equality operator
      EIG_PQ &operator=(EIG_PQ &);

      double *operator[](int);

      double min();

      double max();

      double center_dev();

   private:

      //double pionter, will be allocated in constructor and the eigenvalues will be stored in it
      double **eig;

      //!number of particles
      int N;

      //!dimension of sp space
      int M;

      //!dimension of tp space
      int n_tp;

      //!total dimension of the EIG_PQ vector.
      int dim;

};

#endif
