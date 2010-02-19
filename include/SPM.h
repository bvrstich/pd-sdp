#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * Deze klasse SPM stelt een single particle matrix voor. Hij erft van de klasse Matrix en breidt die uit met wat 
 * specifieke memberfuncties en een kennis van deeltjesaantal en aantal orbitalen.
 */

class SPM : public Matrix {

   /**
    * Output stream operator overloaded, het gebruik is simpel: wil je naar een file printen, maak dan een
    * ifstream object aan en doe \n\n
    * object << spm_p << endl;\n\n
    * Wil je naar het scherm printen:\n\n
    * cout << spm_p << endl;\n\n
    * @param output de stream waarnaar je toe schrijft.
    * @param spm_p de te printen matrix
    */
   friend ostream &operator<<(ostream &output,SPM &spm_p);

   public:
      
      //constructor
      SPM(int M,int N);

      //copy constructor
      SPM(SPM &);

      //destructor
      virtual ~SPM();

      using Matrix::operator=;

      int gN();

      int gM();

      /**
       * construeert een SPM uit een TPM of PHM. Definitie van deze bewerkingen staan in andere nota's.\n
       * Eigenlijk is deze functie bar * scale.
       * @param scale de factor waarmee de bar(MT) vermendigvuldigd wordt vb. 1/(N - 1) in het geval van een 
       * gewone eendeeltjesdichtheidsmatrix.
       * @param MT de PHM of TPM inputmatrix.
       */
      template<class MatrixType>
         void constr(double scale,MatrixType &MT){

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

                  (*this)(a,b) *= scale;

               }

            this->symmetrize();

         }

      /**
       * Constructor van een SPM object met initialisatie door middel van de functie constr .
       * @param scale de factor waarmee de bar(MT) vermendigvuldigd wordt vb. 1/(N - 1) in het geval van een 
       * gewone eendeeltjesdichtheidsmatrix.
       * @param MT de PHM of TPM inputmatrix.
       */
      template<class MatrixType>
         SPM(double scale,MatrixType &MT) : Matrix(MT.gM()) {

            this->M = MT.gM();
            this->N = MT.gN();

            this->constr(scale,MT);

         }

      /**
       * construeert een SPM uit een TPM of PHM. Definitie van deze bewerkingen staan in andere nota's.\n
       * @param MT de PHM of TPM inputmatrix.
       */
      template<class MatrixType>
         void bar(MatrixType &MT){

            for(int a = 0;a < M;++a)
               for(int b = a;b < M;++b){

                  (*this)(a,b) = 0.0;

                  for(int l = 0;l < M;++l)
                     (*this)(a,b) += MT(a,l,b,l);

               }

            this->symmetrize();

         }

   private:

      //!dimensie van de eendeeltjesruimte, tevens ook dimensie van de matrix
      int M;

      //!Aantal deeltjes
      int N;

};

#endif
