#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 09-03-2010\n\n
 * This class, EIG is a "block"-vector over the carrierspace's of the active condtions. It contains room
 * to store the eigenvalues and special member function that work with these eigenvalues.
 * This class should only be used when a SUP matrix has been diagonalized, some functions could give strange results when the EIG object is filled
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
   friend ostream &operator<<(ostream &output,EIG &eig_p);

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


      double operator()(int,int);

      //overload equality operator
      EIG &operator=(EIG &);

      double *operator[](int);

      double min();

      double max();

      double center_dev();

#ifdef __G_CON

      int gn_ph();

#endif

   private:

      //!single pointer to doubles, the eigenvalues of the SUP matrix will be stored here.
      double *eig;

      //!number of particles
      int N;

      //!dimension of sp space
      int M;

      //!dimension of tp space
      int n_tp;

#ifdef __G_CON
      
      //!dimension of ph space
      int n_ph;

#endif

      //!total dimensie of the EIG object
      int dim;

};

#endif
