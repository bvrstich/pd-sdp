#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class SPM : public Matrix {

   friend ostream &operator<<(ostream &,SPM &);

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

      //definitie van template functies mag niet buiten de klassedefinitie
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

      //tconstructor met initialisatie op trace van een template matrix:
      template<class MatrixType>
         SPM(double scale,MatrixType &MT) : Matrix(MT.gM()) {

            this->M = MT.gM();
            this->N = MT.gN();

            this->constr(scale,MT);

         }

      //definitie van template functies mag niet buiten de klassedefinitie
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

      int M;//dim sp space
      int N;//nr of particles

};

#endif
