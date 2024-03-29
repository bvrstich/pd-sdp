#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "include.h"

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

      //default constructor
      EIG();
   
      //constructor met initialisatie op 
      EIG(const SUP &);
      
      //copy constructor
      EIG(const EIG &);

      //destructor
      ~EIG();

      void diagonalize(const SUP &);

      int gN() const;

      int gM() const;

      int gdim() const;

      double centerpot(double,const EIG &,double,double) const;

      //overload equality operator
      EIG &operator=(const EIG &);

      Vector<TPM> &tpv(int) const;

#ifdef __G_CON

      Vector<PHM> &phv() const;

#endif

#ifdef __T1_CON

      Vector<DPM> &dpv() const;

#endif

#ifdef __T2_CON

      Vector<PPHM> &pphv() const;

#endif

#ifdef __T2P_CON

      Vector<T2PM> &t2pv() const;

#endif

      double min() const;

      double max() const;

      double center_dev() const;

      static void init(int,int);

   private:

      //!variable that tells if the memory has been allocated (flag = 1) or not (flag = 0)
      int flag;

      //!double pointer to a Vector<TPM> object, the eigenvalues of the P and Q part of a SUP matrix will be stored here.
      Vector<TPM> **v_tp;

      //!number of particles
      static int N;

      //!dimension of sp space
      static int M;

#ifdef __G_CON

      //!pointer to a Vector<PHM> object that will contain the eigenvalues of the G part of a SUP matrix
      Vector<PHM> *v_ph;
      
#endif

#ifdef __T1_CON

      //!pointer to a Vector<DPM> object that will contain the eigenvalues of the T1 part of a SUP matrix
      Vector<DPM> *v_dp;

#endif

#ifdef __T2_CON

      //!pointer to a Vector<PPHM> object that will contain the eigenvalues of the T2 part of a SUP matrix
      Vector<PPHM> *v_pph;

#endif

#ifdef __T2P_CON

      //!pointer to a Vector<T2PM> object that will contain the eigenvalues of the T2P part of a SUP matrix
      Vector<T2PM> *v_t2p;

#endif

      //!total dimension of the EIG object
      static int dim;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
