#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"

#ifndef PQ

#include "PHM.h"

#endif

class EIG;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * De klasse SUP (supermatrix) is een blokmatrix over de carrierspace's van de verschillende condities:\n
 * Belangrijk om hierbij te onthouden dat dit een algemene blokmatrix is met twee TPM's en een PHM die niet 
 * noodzakelijk met elkaar in verband staan. Zie ook notes primal_dual.pdf .
 */

class SUP{
  
   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << sup_p << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << sup_p << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param sup_p de te printen SUP-matrix
    */
   friend ostream &operator<<(ostream &output,SUP &sup_p);

   public:

      //constructor
      SUP(int M,int N);

      //copy constructor
      SUP(SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(SUP &);

      //overload -= operator
      SUP &operator-=(SUP &);

      //overload equality operator
      SUP &operator=(SUP &);

      //overload equality operator
      SUP &operator=(double &);

      TPM &tpm(int i);

#ifndef PQ

      PHM &phm();

#endif

      //initialiseer S
      void init_S();

      //initialiseer Z
      void init_Z(double alpha,TPM &ham,SUP &u_0);

      int gN();

      int gM();

      int gn_tp();

#ifndef PQ
      
      int gn_ph();

#endif
      int gdim();

      double ddot(SUP &);

      void invert();

      void dscal(double alpha);

      void proj_U();

      //maak de matrix D, nodig voor de hessiaan van het stelsel
      void D(SUP &S,SUP &Z);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP &,SUP &);

      void daxpy(double alpha,SUP &);

      double trace();

      double U_trace();

      void proj_C();

      SUP &mprod(SUP &,SUP &);

      void fill(TPM &);

      int solve(SUP &B,SUP &D);

      void H(SUP &B,SUP &D);

      void proj_U_Tr();

      void diagonalize(EIG &);

      double center_dev(SUP &Z);

      double line_search(SUP &DZ,SUP &S,SUP &Z,double max_dev);

   private:

      //!dubbele pointer van TPM's, deze zullen verwijzen naar de twee TPM matrices
      TPM **SZ_tp;

      //!aantal sp orbitalen
      int M;

      //!aantal deeltjes
      int N;

      //!diminsie van tp space
      int n_tp;

      //!totale dimensie van de supermatrix
      int dim;

#ifndef PQ
      
      //!pointer naar de particle hole matrix
      PHM *SZ_ph;

      //!dimensie van de particle hole ruimte
      int n_ph;

#endif
      
};

#endif
