#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class SPM;
class SUP;
class PHM;
class DPM;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * Deze klasse, TPM, stelt een two particle matrix voor. Hij erft van de klasse Matrix en breidt die uit met een
 * lijst van hoe de tweedeeltjesbasis opgebouwd is uit de eendeeltjesbasis. Er zijn ook 
 * specifieke memberfuncties toegevoegd.
 */

class TPM : public Matrix {

   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << tpm_p << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << tpm_p << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param tpm_p de te printen matrix
    */
   friend ostream &operator<<(ostream &output,TPM &tpm_p);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(TPM &);

      //destructor
      virtual ~TPM();

      using Matrix::operator=;

      using Matrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int a,int b,int c,int d) const;

      //geef N terug
      int gN();

      //geef N terug
      int gM();

      //geef n terug
      int gn();

      void hubbard(double U);

      void collaps(SUP &B);

      //Q afbeelding en zijn inverse
      void Q(int option,TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(TPM &b,SUP &D);

      //los het stelsel op
      int solve(TPM &b,SUP &D);

      //de G down en inverse G up
      void G(int option,PHM &);

      void bar(DPM &);

      void T(int option,DPM &);

      void min_unit(double scale);

      void min_qunit(double scale);

   private:


      //!static lijst, dubbele pointer (n_tp*2) van integers die een tp index i neemt en dan twee sp indices teruggeeft: t2s(i,0) = a  t2s(i,1) = b
      static int **t2s;

      //! static lijst, dubbele pointer (M*M) van integers die twee sp indices neem (a,b) een een tp index teruggeeft (i)
      static int **s2t;

      //!static counter, telt het aantal TPM objecten momenteel aangemaakt in het programma omdat het geheugen waar de t2s en s2t pointers naar wijzen maar 1 maal aangemaakt wordt (static dus)
      static int counter;

      //!aantal deeltjes
      int N;

      //!aantal sp orbital, dimensie van eendeeltjesruimte
      int M;

      //!dimensie van de tweedeeltjesruimte
      int n;

};

#endif
