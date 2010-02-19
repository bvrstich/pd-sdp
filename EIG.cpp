#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "EIG.h"
#include "SUP.h"
#include "lapack.h"

/**
 * standard constructor\n
 * Alloceert twee array's met dimensie n_tp en een met dimensie n_ph en zorg ervoor dat de daartoe voorziene pointers ernaar verwijzen.
 * @param M aantal sp orbitals
 * @param N aantal deeltjes
 */
EIG::EIG(int M,int N){
   
   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

}

/**
 * Copy constructor\n
 * Alloceert twee array's met dimensie n_tp en een met dimensie n_ph en zorg ervoor dat de daartoe voorziene pointers ernaar verwijzen. 
 * De inhoud van eig_c wordt dan gekopieerd naar this.
 * @param eig_c De EIG waarvan de inhoud gekopieerd zal worden naar this
 */
EIG::EIG(EIG &eig_c){

   this->N = N;
   this->M = M;
   this->n_tp = M*(M - 1)/2;

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

}
 
/**
 * overload equality operator
 * @param eig_c Deze EIG zal gekopieerd worden in this.
 */
EIG &EIG::operator=(EIG &eig_c){

   int inc = 1;

   dcopy_(&dim,eig_c.eig[0],&inc,eig[0],&inc);

   return *this;

}

/**
 * Toegang tot de afzonderlijke vectoren.
 * @param i index die aangeeft welke vector je terugkrijgt: i = 0: vector van het bovenste blok(TPM), i = 1 vector van het middenste blok(TPM), i = 2 vector van het onderste blok (PHM).
 * @return de overeenkomstige pointer
 */
double *EIG::operator[](int i){

   return eig[i];

}

/**
 * Constructor op basis van een SUP matrix. Alloceert twee array's met dimensie n_tp en een met dimensie n_ph en zorg ervoor dat de daartoe voorziene pointers ernaar verwijzen. 
 * This wordt dan opgevult met de eigenwaarden van de SUP matrix SZ, opgelet, SZ wordt hierbij vernield en bevat zijn eigenvectoren in de kolommen.
 * @param SZ de bewuste SUP matrix
 */
EIG::EIG(SUP &SZ){

   this->N = SZ.gN();
   this->M = SZ.gM();
   this->n_tp = SZ.gn_tp();

#ifndef PQ

   this->n_ph = M*M;

   dim = 2*n_tp + n_ph;

   eig = new double * [3];
   eig[0] = new double [dim];

   for(int i = 1;i < 3;++i)
      eig[i] = eig[i - 1] + n_tp;

#else
   
   dim = 2*n_tp;

   eig = new double * [2];
   eig[0] = new double [dim];

   for(int i = 1;i < 2;++i)
      eig[i] = eig[i - 1] + n_tp;

#endif

   //dit vernietigd de originele matrix!
   for(int i = 0;i < 2;++i)
      (SZ.tpm(i)).diagonalize(eig[i]);

#ifndef PQ

   (SZ.phm()).diagonalize(eig[2]);

#endif

}

/**
 * Destructor, deallocatie van het geheugen.
 */
EIG::~EIG(){

   delete [] eig[0];
   delete [] eig;

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[0][i] << std::endl;

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_tp;++i)
      std::cout << i << "\t" << eig_p.eig[1][i] << std::endl;

#ifndef PQ

   std::cout << std::endl;

   for(int i = 0;i < eig_p.n_ph;++i)
      std::cout << i << "\t" << eig_p.eig[2][i] << std::endl;

#endif

   return output;

}

/**
 * operator () overloaded\n
 * Toegang van buitenaf tot de elementen in de array
 * @param block in welk blok zit het element dat je wil
 * @param index op welke plaats in dat blok staat het element dat je wil
 * @return het gewenste getal eig[block][index]
 */
double EIG::operator()(int block,int index){

   return eig[block][index];

}

/**
 * @return het minimum van de elementen in de EIG\n
 * Opgelet, dit werkt enkel als de EIG gevuld is door diagonalisatie van een SUP.
 */
double EIG::min(){

   double ward = eig[0][0];

   if(ward > eig[1][0])
      ward = eig[1][0];

#ifndef PQ

   if(ward > eig[2][0])
      ward = eig[2][0];

#endif

   return ward;

}

/**
 * @return het maximum van de elementen in de EIG\n
 * Opgelet, dit werkt enkel als de EIG gevuld is door diagonalisatie van een SUP.
 */
double EIG::max(){

   double ward = eig[0][n_tp - 1];

   if(ward < eig[1][n_tp - 1])
      ward = eig[1][n_tp - 1];

#ifndef PQ

   if(ward < eig[2][n_ph - 1])
      ward = eig[2][n_ph - 1];

#endif

   return ward;

}

/**
 * @return de afwijking van het centraal pad berekend met de logaritmische potentiaalbarriere (zie notes)
 */
double EIG::center_dev(){

   double sum = 0.0;

   for(int i = 0;i < n_tp;++i)
      sum += eig[0][i];

   for(int i = 0;i < n_tp;++i)
      sum += eig[1][i];

   double log_product = 0.0;

   for(int i = 0;i < n_tp;++i)
      log_product += log(eig[0][i]);

   for(int i = 0;i < n_tp;++i)
      log_product += log(eig[1][i]);

#ifndef PQ

   for(int i = 0;i < n_ph;++i)
      sum += eig[2][i];

   for(int i = 0;i < n_ph;++i)
      log_product += log(eig[2][i]);

#endif

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * Geeft terug wat de afwijking is van het centraal pad (grootte van de potentiaal) wanneer je stapgrootte 
 * alpha zet langs de primal dual Newton richting, d.m.v. de eigenwaarden eigen_S en eigen_Z berekend in SUP::line_search.\n
 * (*this) = eigen_S --> de eigenwaarden voor de DS stap
 * @param alpha afstand langs de Newton richting
 * @param eigen_Z  --> de eigenwaarden voor de DZ stap
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter berekend in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter berekend in SUP::line_search
 * @return de afwijking van het centraal pad bij stapgrootte alpha langs de Newton richting
 */
double EIG::centerpot(double alpha,EIG &eigen_Z,double c_S,double c_Z){

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eig[0][i]);

   for(int i = 0;i < dim;++i)
      ward -= log(1.0 + alpha*eigen_Z.eig[0][i]);

   return ward;

}
