#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ofstream;

#include "lapack.h"
#include "Matrix.h"
#include "SPM.h"
#include "TPM.h"

#ifndef PQ

#include "PHM.h"

#endif 

#include "SUP.h"
#include "EIG.h"

int main(void){

   cout.precision(10);

   int M = 12;//dim sp hilbert space
   int N = 5;//nr of particles

   int n_tp = M*(M - 1)/2;//dim van tp ruimte

   int dim = 2*n_tp;

#ifndef PQ

   int n_ph = M*M;//dim van ph ruimte
   dim += n_ph;

#endif

   //hamiltoniaan
   TPM ham(M,N);
   ham.hubbard(1.0);

   SUP S(M,N);
   S.init_S();

   SUP Z(M,N);
   Z.init_Z(10.0,ham,S);

   //eerste primal dual gap:
   double pd_gap = S.ddot(Z);
   double energy = (S.tpm(0)).ddot(ham);

   double center_dev = S.center_dev(Z);

   //eerst centering
   double gamma = 1.0;

   double tolerance = 1.0e-4;

   //flag == 0 : initiele centering run (tot op tolerance)
   //flag == 1 : doe een stap met gamma = 0
   //flag == 2 : doe een stap met gamma = 1
   //flag == 3 : game over man

   int flag = 0;

   double a;//stapgrootte

   while(flag != 3){

      cout << (S.tpm(0)).trace() << "\t" << pd_gap << "\t" << center_dev << "\t" << energy << "\t";

      //matrix D aanmaken voor de hessiaan van het duale stelsel
      SUP D(M,N);
      D.D(S,Z);

      //D inverteren voor de hessiaan van het primale stelsel
      SUP D_inv(D);
      D_inv.invert();

      //rechterlid maken van stelsel dat moet worden opgelost:
      SUP B(S);

      //invert B
      B.invert();

      //schalen met 
      B.dscal(gamma*pd_gap/dim);

      B -= Z;

      //nu kan het rechterlid worden gemaakt:
      TPM b(M,N);

      b.Q(1,B.tpm(1));

      b += B.tpm(0);

#ifndef PQ

      TPM hulp(M,N);
      hulp.G(1,B.phm());

      b += hulp;

#endif

      b.proj_Tr();

      //dit wordt de stap:
      TPM delta(M,N);

      //los het stelsel op, geeft aantal iteraties nodig terug:
      cout << delta.solve(b,D_inv) << "\t";

      //nog updaten van S en Z
      SUP DS(M,N);

      DS.fill(delta);

      //DZ is B - D^{-1}*DS*D^{-1}
      SUP DZ(B);

      //eerst D^{-1}*DS*D^{-1} in DZ stoppen
      B.L_map(D_inv,DS);

      DZ -= B;

      //voor de zekerheid nog projecteren op juiste subruimte:
      DZ.proj_C();

      //met deze 'ansatz' het Z stelsel proberen op te lossen
      //eerst rechterlid B maken
      B = Z;

      B.invert();

      B.dscal(gamma*pd_gap/dim);

      B -= S;

      B.proj_C();

      //los het stelsel op, geeft aantal duale iteraties nodig terug:
      cout << DZ.solve(B,D) << endl;

      //welke stapgrootte moet ik nemen?

      if(flag == 0 || flag == 2){//voor centering

         S += DS;
         Z += DZ;

      }
      else{

         //zoek de ideale afstand (geef ook een waarde mee voor de maximale afwijking van het centraal pad):
         a = DS.line_search(DZ,S,Z,5.0);

         S.daxpy(a,DS);
         Z.daxpy(a,DZ);

      }

      //update van enkele belangrijke variabelen
      pd_gap = S.ddot(Z);
      energy = (S.tpm(0)).ddot(ham);
      center_dev = S.center_dev(Z);

      //keuze voor volgende iteratie:
      if(flag == 0){

         //als hij voldoende gecenterd is, exit.
         if(center_dev < tolerance){

            flag = 1;
            gamma = 0.0;

         }

      }
      else if(flag == 1){

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//center when not convergence

            flag = 2;
            gamma = 1.0;

         }

      }
      else{//flag == 2: dus na een centering stap

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//take another step downwards when not converged

            flag = 1;
            gamma = 0;

         }

      }

   }

   cout << endl;
   cout << "FINAL RESULT " << endl;
   cout << endl;
   cout << "E_0 = " << energy << " with accuracy of " << pd_gap << " and a deviation from centrality of " << center_dev << endl;
   cout << endl;

   return 0;

}
