/**
 * @file 
 * This is a wrapper class around the different SUP_PQ(GT1) classes. It is decided at compile time from which 
 * SUP_* file this class inherits. Compile with PQ to inherit from SUP_PQ, compile with PQG to inherit from SUP_PQG, etc. .
 * This way, you can use the SUP object everywhere in the program without having to worry about which conditions are used.
 */
#ifndef SUP_H
#define SUP_H

//if PQ is defined, inherit from SUP_PQ
#ifdef PQ

#include "SUP/SUP_PQ.h"

class SUP : public SUP_PQ{

   public :

      SUP(int M,int N) : SUP_PQ(M,N) { }

      SUP(SUP &SZ) : SUP_PQ(SZ) { }

      ~SUP(){ }

};

#endif

#endif
