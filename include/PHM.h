#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "Matrix.h"
#include "TPM.h"
#include "SPM.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * Deze klasse, PHM, stelt een particle hole matrix voor. Hij erft van de klasse Matrix en breidt die uit met een
 * lijst van hoe de deeltje-gat-basis opgebouwd is uit de eendeeltjesbasis. Er zijn ook specifieke memberfuncties toegevoegd.
 */

class PHM : public Matrix {

   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << phm_p << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << phm_p << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param phm_p de te printen matrix
    */
   friend ostream &operator<<(ostream &output,PHM &phm_p);

   public:
      
      //constructor
      PHM(int M,int N);

      //copy constructor
      PHM(PHM &);

      //destructor
      virtual ~PHM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      //double operator()(int a,int b,int c,int d) const;

      //change the numbers in sp mode
      double &operator()(int a,int b,int c,int d);

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      //geef dim terug
      int gn();

      void G(int option,TPM &);

      double skew_trace();

      void min_gunit(double scale);

   private:

      //!static counter, telt het aantal PHM objecten momenteel aangemaakt in het programma omdat het geheugen waar de ph2s en s2ph pointers naar wijzen maar 1 maal aangemaakt wordt (static dus)
      static int counter;

      //!static lijst, dubbele pointer (n_ph*2) van integers die een ph index i neemt en dan twee sp indices teruggeeft: ph2s(i,0) = a  ph2s(i,1) = b
      static int **ph2s;

      //! static lijst, dubbele pointer (M*M) van integers die twee sp indices neemt (a,b) een een ph index teruggeeft (i)
      static int **s2ph;

      //!aantal deeltjes
      int N;

      //!aantal sp orbital, dimensie van eendeeltjesruimte
      int M;

      //!dimensie van de particle hole ruimte
      int n;

};

#endif
