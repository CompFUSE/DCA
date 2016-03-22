
//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for LAPACK                                                   *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Peter Staar, ETHZ, December 2011                                         *
// ****************************************************************************

#ifndef LAPACK_C_WRAPPERS
#define LAPACK_C_WRAPPERS

/** \file LAPACK
 *  \author Thomas C. Schulthess, Michael S. Summers, Peter Staar
 */

namespace LAPACK {

  typedef int ftnlen;

// ============================================================================

  extern "C" int  ilaenv_(int*,char*,char*,int*,int*,int*,int*,int nameLen,int optLen);

  extern "C" float  slamch_( char* JOBZ);
  extern "C" double dlamch_( char* JOBZ);

  extern "C" void slaset_(char* UPLO, int* M, int* N, float*  ALPHA, float*  BETA, float*  A, int* LDA);
  extern "C" void dlaset_(char* UPLO, int* M, int* N, double* ALPHA, double* BETA, double* A, int* LDA);
  
  
  // solve

  extern "C" void sgesv_(int*,int*,float*,int*,int*,float*,int*,int*);
  extern "C" void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
  extern "C" void cgesv_(int*,int*,std::complex<float>*,int*,int*,std::complex<float>*,int*,int*);
  extern "C" void zgesv_(int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,int*,int*);

  extern "C" void sgetrs_( char* TRANS,  int *N,  int *NRHS, float *A,  int *LDA, int *IPIV, float *B,  int *LDB, int *INFO);
  extern "C" void dgetrs_( char* TRANS,  int *N,  int *NRHS, double *A,  int *LDA, int *IPIV, double *B,  int *LDB, int *INFO);
  extern "C" void cgetrs_( char* TRANS,  int *N,  int *NRHS, std::complex<float> *A,  int *LDA, int *IPIV, std::complex<float> *B,  int *LDB, int *INFO);
  extern "C" void zgetrs_( char* TRANS,  int *N,  int *NRHS, std::complex<double> *A,  int *LDA, int *IPIV, std::complex<double> *B,  int *LDB, int *INFO);

  // LU factorization

  extern "C" int  sgetrf_(int* M, int* N, float* A, int* LDA, int* IPIV, int* INFO);
  extern "C" int  dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
  extern "C" int  cgetrf_(int* M, int* N, std::complex<float>* A, int* LDA, int* IPIV, int* INFO);
  extern "C" int  zgetrf_(int* M, int* N, std::complex<double>* A, int* LDA, int* IPIV, int* INFO);

  extern "C" int  sgetri_(int* N, float* A, int* LDA, int* IPIV,  float* WORK, int* LWORK, int* INFO);
  extern "C" int  dgetri_(int* N, double* A, int* LDA, int* IPIV,  double* WORK, int* LWORK, int* INFO);
  extern "C" int  cgetri_(int* N, std::complex<float>* A, int* LDA, int* IPIV,  std::complex<float>* WORK, int* LWORK, int* INFO);
  extern "C" int  zgetri_(int* N, std::complex<double>* A, int* LDA, int* IPIV,  std::complex<double>* WORK, int* LWORK, int* INFO);



  extern "C" int  dgeco_(double* A, int* LDA, int* N, int* IPIV, double* rcond, double* z);
  extern "C" int  dgefa_(double* A, int* LDA, int* N, int* IPIV, int* INFO);
  extern "C" int  dgesl_(double* A, int* LDA, int* N, int* IPIV, double* B, int* job);

  extern "C" int  sgeco_(float* A,  int* LDA, int* N, int* IPIV, float* rcond, float* z);
  extern "C" int  sgefa_(float* A,  int* LDA, int* N, int* IPIV, int* INFO);
  extern "C" int  sgesl_(float* A,  int* LDA, int* N, int* IPIV, float* B, int* job);

  // eigenvalue decomposition
  
  /**************************
   ***   GEEV-routines
   **************************/

  extern "C" void sgeev_(char* JOBVL,  char* JOBVR, int* N, float* A,  int* LDA, float* WR, float* WI, float* VL, int* LDVL, float* VR, int* LDVR, float* WORK, int* LWORK,int*   INFO );
  extern "C" void dgeev_(char* JOBVL,  char* JOBVR, int* N, double* A,  int* LDA, double* WR, double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int*   INFO );
  extern "C" void cgeev_(char* JOBVL,  char* JOBVR, int* N, std::complex<float>* A,  int* LDA, std::complex<float>* W, std::complex<float>* VL, int* LDVL, std::complex<float>* VR, int* LDVR, std::complex<float>* WORK, int* LWORK, std::complex<float>* RWORK, int*   INFO );
  extern "C" void zgeev_(char* JOBVL,  char* JOBVR, int* N, std::complex<double>* A,  int* LDA, std::complex<double>* W, std::complex<double>* VL, int* LDVL, std::complex<double>* VR, int* LDVR, std::complex<double>* WORK, int* LWORK, std::complex<double>* RWORK, int*   INFO );

  /**************************
   ***   HEEV-routines
   **************************/

  extern "C" void ssyev_( char* JOBZ,  char* UPLO,  int* N, float* Matrix,  int* LDA, float* eigenvalues, float* WORK,  int* LWORK, int* INFO );
  extern "C" void dsyev_( char* JOBZ,  char* UPLO,  int* N, double* Matrix,  int* LDA, double* eigenvalues, double* WORK,  int* LWORK, int* INFO );
  extern "C" void cheev_( char* jobz,  char* uplo,  int* n, std::complex<float>* a,  int* lda, float* w, std::complex<float>* work,  int* lwork, float* rwork, int* info );
  extern "C" void zheev_( char* jobz,  char* uplo,  int* n, std::complex<double>* a,  int* lda, double* w, std::complex<double>* work,  int* lwork, double* rwork, int* info );

  /**************************
   ***   HEEVD-routines
   **************************/

  extern "C" void ssyevd_( char* JOBZ,  char* UPLO,  int* N, float*  Matrix,  int* LDA, float* eigenvalues , float* WORK , int* LWORK, int* IWORK , int* LIWORK, int* INFO );
  extern "C" void dsyevd_( char* JOBZ,  char* UPLO,  int* N, double* Matrix,  int* LDA, double* eigenvalues, double* WORK, int* LWORK, int* IWORK , int* LIWORK, int* INFO );

  extern "C" void cheevd_( char* jobz,  char* uplo,  int* n, std::complex<float>*  a,  int* lda, float*  w, std::complex<float>*  work,  int* lwork, float*  rwork, int* LRWORK, int* IWORK , int* LIWORK, int* info );
  extern "C" void zheevd_( char* jobz,  char* uplo,  int* n, std::complex<double>* a,  int* lda, double* w, std::complex<double>* work,  int* lwork, double* rwork, int* LRWORK, int* IWORK , int* LIWORK, int* info );

  /**************************
   ***   HEEVX-routines
   **************************/

  extern "C" void ssyevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, float * A, int* LDA, 
			   float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M,
			   float* W, float* Z, int* LDZ, 
			   float* WORK,  int* LWORK, int* IWORK, int* IFAIL, int* INFO );

  extern "C" void dsyevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, double* A, int* LDA, 
			   double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
			   double* W, double* Z, int* LDZ, 
			   double* WORK,  int* LWORK, int* IWORK, int* IFAIL, int* INFO );

  extern "C" void cheevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, std::complex<float>* A, int* LDA, 
			   float* VL, float* VU, int* IL, int* UL, 
			   float* ABSTOL, int* M,
			   float* W, std::complex<float>* Z, int* LDZ, 
			   std::complex<float>* WORK,  int* LWORK, float* RWORK,  int* LRWORK, int* IWORK, int* IFAIL, int* INFO);

  extern "C" void zheevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, std::complex<double>* A, int* LDA, 
			   double* VL, double* VU, int* IL, int* UL, 
			   double* ABSTOL, int* M,
			   double* W, std::complex<double>* Z, int* LDZ, 
			   std::complex<double>* WORK,  int* LWORK, double* RWORK,  int* LRWORK, int* IWORK, int* IFAIL, int* INFO);

  /**************************
   ***   HEEVR-routines
   **************************/

  extern "C" void ssyevr_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, float * A, int* LDA, 
			   float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M,
			   float* W, float* Z, int* LDZ, int* ISUPPZ, 
			   float* WORK,  int* LWORK, int* IWORK, int* LIWORK, int* INFO );

  extern "C" void dsyevr_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, double* A, int* LDA, 
			   double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
			   double* W, double* Z, int* LDZ, int* ISUPPZ, 
			   double* WORK,  int* LWORK, int* IWORK, int* LIWORK, int* INFO );

  extern "C" void cheevr_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, std::complex<float>* A, int* LDA, 
			   float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M,
			   float* W, std::complex<float>* Z, int* LDZ, int* ISUPPZ, 
			   std::complex<float>* WORK, int* LWORK, float* RWORK,  int* LRWORK, int* IWORK, int* LIWORK, int* INFO );

  extern "C" void zheevr_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, std::complex<double>* A, int* LDA, 
			   double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
			   double* W, std::complex<double>* Z, int* LDZ, int* ISUPPZ, 
			   std::complex<double>* WORK,  int* LWORK, double* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO );

  /**************************
   ***   GESVD-routines
   **************************/
  
  extern "C" void sgesvd_(char* JOBU, char* JOBVT, int* M, int* N, float* A              , int* LDA, float* S, float* U              , int* LDU, float* VT              , int* LDVT, float* WORK, int* LWORK, int* INFO);
  extern "C" void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A              , int* LDA, double* S, double* U              , int* LDU, double* VT              , int* LDVT, double* WORK, int* LWORK, int* INFO);

  extern "C" void cgesvd_(char* JOBU, char* JOBVT, int* M, int* N, std::complex<float>* A, int* LDA, float* S, std::complex<float>* U, int* LDU, std::complex<float>* VT, int* LDVT, std::complex<float>* WORK, int* LWORK, float* RWORK, int* INFO );
  extern "C" void zgesvd_(char* JOBU, char* JOBVT, int* M, int* N, std::complex<double>* A, int* LDA, double* S, std::complex<double>* U, int* LDU, std::complex<double>* VT, int* LDVT, std::complex<double>* WORK, int* LWORK, double* RWORK, int* INFO );

  /**************************
   ***   GESVDD-routines
   **************************/

  extern "C" void sgesdd_(char* JOBZ, int* M, int* N, float* A              , int* LDA, float* S, float* U              , int* LDU, float* VT              , int* LDVT, float* WORK, int* LWORK, int* IWORK, int* INFO);
  extern "C" void dgesdd_(char* JOBZ, int* M, int* N, double* A              , int* LDA, double* S, double* U              , int* LDU, double* VT              , int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);

  extern "C" void cgesdd_(char* JOBZ, int* M, int* N, std::complex<float>* A, int* LDA, float* S, std::complex<float>* U, int* LDU, std::complex<float>* VT, int* LDVT, std::complex<float>* WORK, int* LWORK, float* RWORK, int* IWORK, int* INFO );
  extern "C" void zgesdd_(char* JOBZ, int* M, int* N, std::complex<double>* A, int* LDA, double* S, std::complex<double>* U, int* LDU, std::complex<double>* VT, int* LDVT, std::complex<double>* WORK, int* LWORK, double* RWORK, int* IWORK, int* INFO );



  // least square fit

  extern "C" void sgelss_(int* M, int* N, int* NRHS,              float*   A, int* LDA,              float*   B, int* LDB, float*  S, float*  RCOND, int* RANK,              float*   WORK, int* LWORK, int* INFO);
  extern "C" void dgelss_(int* M, int* N, int* NRHS,              double*  A, int* LDA,              double*  B, int* LDB, double* S, double* RCOND, int* RANK,              double*  WORK, int* LWORK, int* INFO);
  extern "C" void cgelss_(int* M, int* N, int* NRHS, std::complex<float>*  A, int* LDA, std::complex<float>*  B, int* LDB, float*  S, float*  RCOND, int* RANK, std::complex<float>*  WORK, int* LWORK, float*  RWORK, int* INFO);
  extern "C" void zgelss_(int* M, int* N, int* NRHS, std::complex<double>* A, int* LDA, std::complex<double>* B, int* LDB, double* S, double* RCOND, int* RANK, std::complex<double>* WORK, int* LWORK, double* RWORK, int* INFO);

  // generalized least square fit
  extern "C" void sgglse_(int* M, int* N, int* P, float*  A, int* LDA, float*  B, int* LDB, float*  C, float*  D, float*  X, float* WORK, int* LWORK, int* INFO);
  extern "C" void dgglse_(int* M, int* N, int* P, double* A, int* LDA, double* B, int* LDB, double* C, double* D, double* X, double* WORK, int* LWORK, int* INFO);

  // QR

  extern "C" void sgeqp3_(int* M, int* N, float* A, int* LDA, int* JPVT, float* TAU, float* WORK, int* LWORK, int* INFO);
  extern "C" void dgeqp3_(int* M, int* N, double* A, int* LDA, int* JPVT, double* TAU, double* WORK, int* LWORK, int* INFO);

  extern "C" void sorgqr_(int* M, int* N, int* K, float*  A, int* LDA, float* TAU,  float* WORK,  int* LWORK, int* INFO);
  extern "C" void dorgqr_(int* M, int* N, int* K, double* A, int* LDA, double* TAU, double* WORK, int* LWORK, int* INFO);

// ============================================================================

/*
  inline int ILAENV(int ispec,
		    std::string name, 
		    std::string opt,
		    int n1, int n2, int n3, int n4) {
    ftnlen optSize = opt.size();
    if (1>optSize) optSize =1;
    return ilaenv_(&ispec,const_cast<char*>(name.c_str()),
		   const_cast<char*>(opt.c_str()),
		   &n1,&n2,&n3,&n4,
		   (ftnlen)name.size(),optSize);
  }
  inline void GESV(int ma,int mb,float* a,int lda,int* pivot,
		   float* b,int ldb,int& info) {
    sgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,double* a,int lda,int* pivot,
		   double* b,int ldb,int& info) {
    dgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,std::complex<float>* a,int lda,int* pivot,
		   std::complex<float>* b,int ldb,int& info) {
    cgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,std::complex<double>* a,int lda,int* pivot,
		   std::complex<double>* b,int ldb,int& info) {
    zgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }

  inline void GETRF(int ma, int na, double* a, int lda, int* pivot, int& info) {
    dgetrf_(&ma,&na,a,&lda,pivot,&info);
  }
  
  inline void GETRF(int ma, int na, float* a, int lda, int* pivot, int& info) {
    sgetrf_(&ma,&na,a,&lda,pivot,&info);
  }

  //======================================================================
  
  inline void GETRI(int na, float* a, int lda, int* pivot, float* work, int lwork, int& info) {
    sgetri_(&na,a,&lda,pivot,work,&lwork,&info);
  }

  inline void GETRI(int na, double* a, int lda, int* pivot, double* work, int lwork, int& info) {
    dgetri_(&na,a,&lda,pivot,work,&lwork,&info);
  }
  
  template<typename T> inline int GETRI_BlockSize(int N);
  template<> inline int GETRI_BlockSize<float>(int N) {
    enum { BLOCKSIZE=1};
    return ILAENV(BLOCKSIZE,"SGETRI","1234",N,-1,-1,-1);

  }
  template<> inline int GETRI_BlockSize<double>(int N) {
    enum { BLOCKSIZE=1};
    return ILAENV(BLOCKSIZE,"DGETRI","1234",N,-1,-1,-1);
  }

  //======================================================================
  
  inline void GECO(double* a, int lda, int n, int* pivot, double& rcond, double* z) {
    dgeco_(a,&lda,&n,pivot,&rcond,z);
  }

  inline void GEFA(double* a, int lda, int n, int* pivot, int& info) {
    dgefa_(a,&lda,&n,pivot,&info);
  }
  
  inline void GESL(double* a, int lda, int n, int* pivot, double* b, int job) {
    dgesl_(a,&lda,&n,pivot,b,&job);
  }
  
  //======================================================================
  
  inline void GECO(float* a, int lda, int n, int* pivot, float& rcond, float* z) {
    sgeco_(a,&lda,&n,pivot,&rcond,z);
  }

  inline void GEFA(float* a, int lda, int n, int* pivot, int& info) {
    sgefa_(a,&lda,&n,pivot,&info);
  }
  
  inline void GESL(float* a, int lda, int n, int* pivot, float* b, int job) {
    sgesl_(a,&lda,&n,pivot,b,&job);
  }
*/
}      /* namespace LAPACK */

#endif // PSIMAG_LAPACK_H
