/**
 * @file 
 * This is a wrapper class around the different EIG_PQ(GT1) classes. It is decided at compile time from which 
 * EIG_* file this class inherits. Compile with PQ to inherit from EIG_PQ, compile with PQG to inherit from EIG_PQG, etc. .
 * This way, you can use the EIG object everywhere in the program without having to worry about which conditions are used.
 */
#ifndef EIG_H
#define EIG_H

//if PQ is defined, inherit from EIG_PQ
#ifdef PQ

#include "EIG/EIG_PQ.h"

class EIG : public EIG_PQ { 

   public :

      EIG(int M,int N) : EIG_PQ(M,N) { }

      EIG(EIG &eig) : EIG_PQ(eig) { }

      EIG(SUP_PQ &SZ) : EIG_PQ(SZ) { }

      ~EIG(){ }
   
};

#endif

#endif
