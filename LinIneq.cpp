#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

int LinIneq::nr;

double LinIneq::a;
double LinIneq::b;
double LinIneq::c;

LinCon **LinIneq::li;

double *LinIneq::coef;

/**
 * Static function that initialializes the static members: it allocates the LinCon objects and 
 * initializes them.
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param nr_in nr of contraint matrices
 */
void LinIneq::init(int M,int N,int nr_in){

   nr = nr_in;

   //allocate
   li = new LinCon * [nr];

   for(int i = 0;i < nr;++i)
      li[i] = new LinCon(M,N);

   //fill
   for(int i = 0;i < nr;++i)
      li[i]->fill_Random();

   //what are the coef's of the overlap matrix without the linear constraints:
   constr_overlap(M,N);

   int n = 2*nr + 1;

   //make the linear system for the inverse overlapmatrix coefficient calculations
   coef = new double [n*n];

   //now make some variables needed for the making of the system

   //traces:
   double con_tr[nr];

   for(int i = 0;i < nr;++i)
      con_tr[i] = li[i]->gI().trace();

   //overlaps
   Matrix I_overlap(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         I_overlap(i,j) = li[i]->gI().ddot(li[j]->gI());

   I_overlap.symmetrize();

   //tenslotte
   Matrix I_bar_overlap(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         I_bar_overlap(i,j) = li[i]->gI_bar().ddot(li[j]->gI_bar());

   I_bar_overlap.symmetrize();

   //fill the system:

   //first column
   coef[0] = a - 2.0*c*(M - 1.0) + b*M*(M - 1.0);

   for(int i = 1;i <= nr;++i)
      coef[i] = 0.0;

   for(int i = nr + 1;i <= 2*nr;++i)
      coef[i] = (b*(M - 1.0) - c)*con_tr[i - nr - 1];

   for(int k = 1;k <= nr;++k){//columns

      coef[k*n] = con_tr[k - 1];

      for(int i = 1;i <= nr;++i){//rows

         coef[k*n + i] = I_overlap(i - 1,k - 1);

         if(i == k)
            coef[k*n + i] += a + b;

      }

      for(int i = nr + 1;i <= 2*nr;++i)//rows
         coef[k*n + i] = I_bar_overlap(i - nr - 1,k - 1);

   }

   for(int k = nr + 1;k <= 2*nr;++k){//columns

      coef[k*n] = 0.0;

      for(int i = 1;i <= nr;++i)//rows
         coef[k*n + i] = 0.0;

      coef[k*n + k - nr] = -4.0*c;

      for(int i = nr + 1;i <= 2*nr;++i)
         coef[k*n + i] = 0.0;

      coef[k*n + k] = a - c*(M - 2.0);

   }

   int *ipiv = new int [n];

   int lwork = 3*n;

   double *work = new double [lwork];

   int info;

   //done, now invert the mofo
   dgetrf_(&n,&n,coef,&n,ipiv,&info);

   dgetri_(&n,coef,&n,ipiv,work,&lwork,&info);

   delete [] ipiv;
   delete [] work;

}

/**
 * delete the statics
 */
void LinIneq::clean(){

   for(int i = 0;i < nr;++i)
      delete li[i];

   delete [] li;

   delete [] coef;

}

/**
 * Will calculate the parameters needed for the overlapmatrix-map: S_a,S_b and S_c.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
void LinIneq::constr_overlap(int M,int N){

   a = 1.0;
   b = 0.0;
   c = 0.0;

#ifdef __Q_CON

   a += 1.0;
   b += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   c += (2.0*N - M)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __G_CON

   a += 4.0;
   c += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   a += M - 4.0;
   b += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   c -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON
   
   a += 5.0*M - 8.0;
   b += 2.0/(N - 1.0);
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

#ifdef __T2P_CON

   a += 5.0*M - 4.0;
   b += 2.0/(N - 1.0);
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M - 2.0)/(2.0*(N - 1.0)*(N - 1.0));

#endif

}

/**
 * constructor of a LinIneq object
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
LinIneq::LinIneq(int M,int N){

   this->M = M;
   this->N = N;

   proj = new double [nr];

   proj_bar = new double [nr];

}

/**
 * copy constructor
 * @param li_copy The LinIneq object to be copied
 */
LinIneq::LinIneq(const LinIneq &li_copy){

   this->M = li_copy.gM();
   this->N = li_copy.gN();

   this->tr = li_copy.gtr();

   proj = new double [nr];

   for(int i = 0;i < nr;++i)
      proj[i] = li_copy.gproj(i);

   proj_bar = new double [nr];

   for(int i = 0;i < nr;++i)
      proj_bar[i] = li_copy.gproj_bar(i);

}

/**
 * destructor
 */
LinIneq::~LinIneq(){

   delete [] proj;
   delete [] proj_bar;

}

/**
 * @return the number of currently applied linear constraints
 */
int LinIneq::gnr() const{

   return nr;

}

/**
 * read and write access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
LinCon &LinIneq::operator[](int i){

   return *li[i];

}

/**
 * read only access to your LinCon object
 * @param i row number
 * @return the entry on index i
 */
const LinCon &LinIneq::operator[](int i) const{

   return *li[i];

}

/**
 * @return nr of particles
 */
int LinIneq::gN() const{

   return N;

}

/**
 * @return nr of sp orbitals
 */
int LinIneq::gM() const{

   return M;

}

/**
 * Calculate the projection of the input TPM object on the constraint matrices of the elements of LinIneq
 * @param tpm The input TPM.
 */
void LinIneq::fill(const TPM &tpm){

   for(int i = 0;i < nr;++i)
      proj[i] = 4.0*(li[i]->gI()).ddot(tpm);

   SPM spm(M,N);
   spm.bar(tpm);

   for(int i = 0;i < nr;++i)
      proj_bar[i] = (li[i]->gI_bar()).ddot(spm);

   tr = 2.0*tpm.trace();

}

/**
 * @param i the index of the constraint we are interested in.
 * @return the projection on the constaint with index i
 */
double LinIneq::gproj(int i) const {

   return proj[i];

}

/**
 * @return the array containing the projection on the constraints.
 */
double *LinIneq::gproj() {

   return proj;

}

/**
 * @param i the index of the constraint we are interested in.
 * @return the barred projection on the barred constaint with index i
 */
double LinIneq::gproj_bar(int i) const {

   return proj_bar[i];

}

/**
 * @return the array containing the barred projection on the barred constraints.
 */
double *LinIneq::gproj_bar(){

   return proj_bar;

}

/**
 * @return the trace of the input TPM scaled with N(N-1)/2
 */
double LinIneq::gtr() const{

   return tr;

}

/**
 * @return the value of the "a" coefficient of the overlapmatrix-map
 */ 
double LinIneq::ga() const {

   return a;

}

/**
 * @return the value of the "b" coefficient of the overlapmatrix-map
 */ 
double LinIneq::gb() const {

   return b;

}

/**
 * @return the value of the "c" coefficient of the overlapmatrix-map
 */ 
double LinIneq::gc() const {

   return c;

}

/**
 * The "alpha" function, see notes, projects a LinIneq (actually a TPM) onto a scalar.
 * @return the alpha function value for this LinIneq object
 */
double LinIneq::alpha() const {

   //I AM HERE!! WRITE ME!!
   return 0.0;

}

/**
 * The "beta" function, see notes, projects a LinIneq (actually a TPM) onto a scalar.
 * @return the beta function value for this LinIneq object
 */
double LinIneq::beta(int index) const {

   //ME TOO!!
   return 0.0;

}
