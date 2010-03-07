#ifndef SUP_PQ_H
#define SUP_PQ_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"

class EIG_PQ;

/**
 * @author Brecht Verstichel
 * @date 04-03-2010\n\n
 * This class, SUP_PQ is a blockmatrix over the carrierspace's of the P and Q conditions. This is the
 * motherclass of all the SUP_PQ* classes because I assume that the P and Q condtional will always be used.
 * This class contains two TPM objects, that are independent of each other (by which I mean that TPM::Q(SUP_PQ::tpm (0))
 * is not neccesarily equal to SUP_PQ::tpm (1)).
 */
class SUP_PQ{
  
    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP_PQ you want to print
    */
   friend ostream &operator<<(ostream &output,SUP_PQ &sup_p);

   public:

      //constructor
      SUP_PQ(int M,int N);

      //copy constructor
      SUP_PQ(SUP_PQ &);

      //destructor
      ~SUP_PQ();

      //overload += operator
      SUP_PQ &operator+=(SUP_PQ &);

      //overload -= operator
      SUP_PQ &operator-=(SUP_PQ &);

      //overload equality operator
      SUP_PQ &operator=(SUP_PQ &);

      //overload equality operator
      SUP_PQ &operator=(double &);

      TPM &tpm(int i);

      //initialiseer S
      void init_S();

      //initialiseer Z
      void init_Z(double alpha,TPM &ham,SUP_PQ &u_0);

      int gN();

      int gM();

      int gn_tp();

      int gdim();

      double ddot(SUP_PQ &);

      void invert();

      void dscal(double alpha);

      void proj_U();

      //maak de matrix D, nodig voor de hessiaan van het stelsel
      void D(SUP_PQ &S,SUP_PQ &Z);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQ &,SUP_PQ &);

      void daxpy(double alpha,SUP_PQ &);

      double trace();

      double U_trace();

      void proj_C();

      SUP_PQ &mprod(SUP_PQ &,SUP_PQ &);

      void fill(TPM &);

      int solve(SUP_PQ &B,SUP_PQ &D);

      void H(SUP_PQ &B,SUP_PQ &D);

      void proj_U_Tr();

      EIG_PQ diagonalize();

      virtual EIG_PQ *get_EIG();

      double center_dev(SUP_PQ &Z);

      double line_search(SUP_PQ &DZ,SUP_PQ &S,SUP_PQ &Z,double max_dev);

   protected:

      //!double pointer of TPM's,
      TPM **SZ_tp;

      //!nr of sp orbitals
      int M;

      //!nr of particles
      int N;

      //!dimension of tp space
      int n_tp;

      //!total dimension of the SUP_PQ matrix
      int dim;

};

#endif
