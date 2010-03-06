#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, EIG_PQG is a "block"-vector over the carrierspace's of the P, Q and G conditions. It inherits from EIG_PQ and expands it with
 * a vector over the G-space (of dimension n_ph). Some functions of EIG_PQ are reimplemented to include the G-condition.
 * This class should only be used when a SUP_PQG matrix has been diagonalized, some functions could give strange results when the EIG object is filled
 * with random numbers.\n\n
 */
class EIG{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << eig_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << eig_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param eig_p the EIG you want to print
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:
      
      //constructor
      EIG(int M,int N);

      //copy constructor
      EIG(EIG &);

      //constructor met initialisatie op 
      EIG(SUP &);

      //destructor
      ~EIG();

      int gN();

      int gM();

      int gn_tp();

      double centerpot(double,EIG &,double,double);

#ifndef PQ

      int gn_ph();

#endif

      double operator()(int,int);

      //overload equality operator
      EIG &operator=(EIG &);

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

#ifndef PQ
      
      //!dimension of ph space
      int n_ph;

#endif

      //!total dimension of the EIG vector.
      int dim;

};

#endif
