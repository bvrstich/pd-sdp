/**
 * @mainpage 
 * This is an implementation of a primal dual interior point method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions.
 * The method used is a path following algorithm with predictor corrector steps.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2, PQGT (T1 and T2), PQGT2P and PQGTP (T1 and T2P).
 * @author Brecht Verstichel, Ward Poelmans
 * @date 16-08-2010
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * 
 * In the main the actual program is run.\n 
 * Part 1: An easy initial point is taken and then centered to the required precision (flag == 0)\n
 * Part 2: When the primal dual point is sufficiently centered steps are taken to reduce the
 * primal dual gap and take a large step in that direction (predictor) (flag == 1)\n
 * After each step a correcting step (flag == 2) is taken that brings the primal dual point closer to
 * the central path.\n
 * Part 3: When the primal dual gap is smaller that the required accuracy exit the while. (flag == 3)\n
 * For more information on the actual method, see primal_dual.pdf
 */
int main(void){

   cout.precision(15);

   int M = 8;//dim sp hilbert space
   int N = 4;//nr of particles

   SPM::init(M,N);
   TPM::init(M,N);

#ifdef __G_CON
   PHM::init(M,N);
#endif

#ifdef __T1_CON
   DPM::init(M,N);
#endif

#ifdef __T2_CON
   PPHM::init(M,N);
#endif

#ifdef __T2P_CON
   T2PM::init(M,N);
#endif
   
   SUP::init(M,N);
   EIG::init(M,N);

   TPM ham;
   ham.hubbard_1D(0,1);

   SUP Z;
   Z.init_Z();

   SUP X;
   X.init_X(1000.0,ham,Z);

   int dim = Z.gdim();

   //eerste primal dual gap:
   double pd_gap = Z.ddot(X);
   double energy = (Z.tpm(0)).ddot(ham);

   double center_dev = Z.center_dev(X);

   //eerst centering
   double gamma = 1.0;

   double tolerance = 1.0e-8;

   //flag == 0 : initiele centering run (tot op tolerance)
   //flag == 1 : doe een stap met gamma = 0
   //flag == 2 : doe een stap met gamma = 1
   //flag == 3 : game over man
   int flag = 0;

   double a;//stapgrootte

   while(flag != 3){

      cout << (Z.tpm(0)).trace() << "\t" << pd_gap << "\t" << center_dev << "\t" << energy << "\t";

      //matrix D aanmaken voor de hessiaan van het duale stelsel
      SUP D;
      D.D(Z,X);

      //D inverteren voor de hessiaan van het primale stelsel
      SUP D_inv(D);
      D_inv.invert();

      //rechterlid maken van stelsel dat moet worden opgelost:
      SUP B(Z);

      //invert B
      B.invert();

      //schalen met 
      B.dscal(gamma*pd_gap/dim);

      B -= X;

      //collaps B onto b to construct the right hand side of the primal Newton equation
      TPM b;

      b.collaps(1,B);

      //dit wordt de stap:
      TPM delta;

      //los het stelsel op, geeft aantal iteraties nodig terug:
      cout << delta.solve(b,D_inv) << "\t";

      //nog updaten van S en Z
      SUP DZ;

      DZ.fill(delta);

      //DZ is B - D^{-1}*DS*D^{-1}
      SUP DX(B);

      //eerst D^{-1}*DS*D^{-1} in DZ stoppen
      B.L_map(D_inv,DZ);

      DX -= B;

      //voor de zekerheid nog projecteren op juiste subruimte:
      DX.proj_C();

      //met deze 'ansatz' het Z stelsel proberen op te lossen
      //eerst rechterlid B maken
      B = X;

      B.invert();

      B.dscal(gamma*pd_gap/dim);

      B -= Z;

      B.proj_C();

      //los het stelsel op, geeft aantal duale iteraties nodig terug:
      cout << DX.solve(B,D) << endl;

      //welke stapgrootte moet ik nemen?
      if(flag == 0 || flag == 2){//voor centering

         Z += DZ;
         X += DX;

      }
      else{

         //zoek de ideale afstand (geef ook een waarde mee voor de maximale afwijking van het centraal pad):
         a = DZ.line_search(DX,Z,X,1.0);

         Z.daxpy(a,DZ);
         X.daxpy(a,DX);

      }

      //update van enkele belangrijke variabelen
      pd_gap = Z.ddot(X);
      energy = (Z.tpm(0)).ddot(ham);
      center_dev = Z.center_dev(X);

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

   //print density matrix to file
//   (S.tpm(0)).out("rdm.out");
   for(int i = 0;i < Z.tpm(0).gn();++i)
      for(int j = i;j < Z.tpm(0).gn();++j)
         cout << i << "\t" << j << "\t" << Z.tpm(0)(i,j) << endl;

#ifdef __T2P_CON
   T2PM::clear();
#endif

#ifdef __T2_CON
   PPHM::clear();
#endif

#ifdef __T1_CON
   DPM::clear();
#endif

#ifdef __G_CON
   PHM::clear();
#endif

   TPM::clear();

   return 0;

}

/* vim: set ts=3 sw=3 expandtab :*/
