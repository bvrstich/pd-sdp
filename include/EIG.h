#ifndef EIG_H
#define EIG_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 19-02-2010\n\n
 * De klasse EIG is een blok-vector over de carrierspace's van de verschillende condities: dus eigenlijk een 
 * klasse om makkelijk met de eigenwaarden van een SUP matrix te kunnen werken.
 * Deze wordt enkel gebruikt bij diagonalisatie van een SUP matrix, moest je de EIG objecten random vullen zouden
 * enkele functies rare(foute) resultaten kunnen geven.
 */

class EIG{

   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << eig_p << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << eig_p << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param eig_p de te printen matrix
    */
   friend ostream &operator<<(ostream &output,const EIG &eig_p);

   public:
      
      //constructor
      EIG(int M,int N);

      //copy constructor
      EIG(EIG &);

      //constructor met initialisatie op 
      EIG(SUP &);

      //destructor
      ~EIG();

      int gN();

      int gM();

      int gn_tp();

      double centerpot(double,EIG &,double,double);

#ifndef PQ

      int gn_ph();

#endif

      double operator()(int,int);

      //overload equality operator
      EIG &operator=(EIG &);

      double *operator[](int);

      double min();

      double max();

      double center_dev();

   private:

      //!dubbele pointer, hier zullen de eigenwaarden van de twee TPM 's en de PHM in opgeslaan worden
      double **eig;

      //!aantal deeltjes
      int N;

      //!dimensie van de sp ruimte
      int M;

      //!dimensie van de tp ruimte
      int n_tp;

#ifndef PQ
      
      //!dimensie van de ph ruimte
      int n_ph;

#endif

      //!totale dimensie van EIG vector
      int dim;

};

#endif
