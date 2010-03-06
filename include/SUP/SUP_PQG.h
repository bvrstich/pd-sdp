#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"
#include "PHM.h"

class EIG;

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, SUP_PQG is a blockmatrix over the carrierspace's of the P, Q and G conditions. It inherits form its motherclass
 * SUP_PQ. This class will expand its mother by a pointer to a PHM object (which is independent of the TPM objects, by which I mean that
 * PHM::G( SUP_PQG::tpm (0)) is not neccesarily equal to SUP_PQG::phm()), and it will also redefine its mothers functions to include the PHM contribution.
 */
class SUP{
  
    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP you want to print
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

      PHM &phm();

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

      //!double pointer of TPM's,
      TPM **SZ_tp;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

      //!dimension of tp space
      int n_tp;

      //!total dimension of the SUP matrix
      int dim;

#ifndef PQ
      
      //!pointer to the PHM
      PHM *SZ_ph;

      //!dimension of ph space.
      int n_ph;

#endif
      
};

#endif
