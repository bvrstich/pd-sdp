#ifndef SUP_PQG_H
#define SUP_PQG_H

#include <iostream> 
#include <fstream> 

using std::ostream;

#include "TPM.h"
#include "PHM.h"

#include "SUP_PQ.h"

class EIG_PQG;

/**
 * @author Brecht Verstichel
 * @date 06-03-2010\n\n
 * This class, SUP_PQG_PQG is a blockmatrix over the carrierspace's of the P, Q and G conditions. It inherits form its motherclass
 * SUP_PQG_PQ. This class will expand its mother by a pointer to a PHM object (which is independent of the TPM objects, by which I mean that
 * PHM::G( SUP_PQG_PQG::tpm (0)) is not neccesarily equal to SUP_PQG_PQG::phm()), and it will also redefine its mothers functions to include the PHM contribution.
 */
class SUP_PQG : public SUP_PQ {
  
    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << sup_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << sup_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param sup_p the SUP_PQG you want to print
    */
   friend ostream &operator<<(ostream &output,SUP_PQG &sup_p);

   public:

      //constructor
      SUP_PQG(int M,int N);

      //copy constructor
      SUP_PQG(SUP_PQG &);

      //destructor
      ~SUP_PQG();

      //overload += operator
      SUP_PQG &operator+=(SUP_PQG &);

      //overload -= operator
      SUP_PQG &operator-=(SUP_PQG &);

      //overload equality operator
      SUP_PQG &operator=(SUP_PQG &);

      //overload equality operator
      SUP_PQG &operator=(double &);

      PHM &phm();

      //initialiseer S
      void init_S();

      //initialiseer Z
      void init_Z(double alpha,TPM &ham,SUP_PQG &u_0);

      int gn_ph();

      double ddot(SUP_PQG &);

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(SUP_PQG &,SUP_PQG &);

      void daxpy(double alpha,SUP_PQG &);

      double trace();

      double U_trace();

      SUP_PQG &mprod(SUP_PQG &,SUP_PQG &);

      void fill(TPM &);

      void proj_U_Tr();

      EIG_PQG diagonalize();

      virtual EIG_PQG *get_EIG();

   private:

      //!pointer to the PHM
      PHM *SZ_ph;

      //!dimension of ph space.
      int n_ph;

};

#endif
