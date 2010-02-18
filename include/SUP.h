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

class SUP{
  
   friend ostream &operator<<(ostream &,SUP &);

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

      TPM **SZ_tp;

      int M;
      int N;
      int n_tp;

      int dim;

#ifndef PQ
      
      PHM *SZ_ph;

      int n_ph;

#endif
      
};

#endif
