//-*-C++-*-

#ifndef LINALG_LAPACK_WRAPPERS
#define LINALG_LAPACK_WRAPPERS

namespace LIN_ALG {

  namespace LAPACK {

    extern "C" void sgetrs_( char* TRANS,  int *N,  int *NRHS, float *A,  int *LDA, int *IPIV, float *B,  int *LDB, int *INFO);
    extern "C" void dgetrs_( char* TRANS,  int *N,  int *NRHS, double *A,  int *LDA, int *IPIV, double *B,  int *LDB, int *INFO);
    extern "C" void cgetrs_( char* TRANS,  int *N,  int *NRHS, std::complex<float> *A,  int *LDA, int *IPIV, std::complex<float> *B,  int *LDB, int *INFO);
    extern "C" void zgetrs_( char* TRANS,  int *N,  int *NRHS, std::complex<double> *A,  int *LDA, int *IPIV, std::complex<double> *B,  int *LDB, int *INFO);
    
    extern "C" int  sgetrf_(int* M, int* N, float* A, int* LDA, int* IPIV, int* INFO);
    extern "C" int  dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
    extern "C" int  cgetrf_(int* M, int* N, std::complex<float>* A, int* LDA, int* IPIV, int* INFO);
    extern "C" int  zgetrf_(int* M, int* N, std::complex<double>* A, int* LDA, int* IPIV, int* INFO);
    
    extern "C" int  sgetri_(int* N, float* A, int* LDA, int* IPIV,  float* WORK, int* LWORK, int* INFO);
    extern "C" int  dgetri_(int* N, double* A, int* LDA, int* IPIV,  double* WORK, int* LWORK, int* INFO);
    extern "C" int  cgetri_(int* N, std::complex<float>* A, int* LDA, int* IPIV,  std::complex<float>* WORK, int* LWORK, int* INFO);
    extern "C" int  zgetri_(int* N, std::complex<double>* A, int* LDA, int* IPIV,  std::complex<double>* WORK, int* LWORK, int* INFO);
    
    // eigenvalue decomposition
    extern "C" void sgeev_(char* JOBVL,  char* JOBVR, int* N, float* A,  int* LDA, float* WR, float* WI, float* VL, int* LDVL, float* VR, int* LDVR, float* WORK, int* LWORK,int*   INFO );
    extern "C" void dgeev_(char* JOBVL,  char* JOBVR, int* N, double* A,  int* LDA, double* WR, double* WI, double* VL, int* LDVL, double* VR, int* LDVR, double* WORK, int* LWORK, int*   INFO );
    extern "C" void cgeev_(char* JOBVL,  char* JOBVR, int* N, std::complex<float>* A,  int* LDA, std::complex<float>* W, std::complex<float>* VL, int* LDVL, std::complex<float>* VR, int* LDVR, std::complex<float>* WORK, int* LWORK, float* RWORK, int*   INFO );
    extern "C" void zgeev_(char* JOBVL,  char* JOBVR, int* N, std::complex<double>* A,  int* LDA, std::complex<double>* W, std::complex<double>* VL, int* LDVL, std::complex<double>* VR, int* LDVR, std::complex<double>* WORK, int* LWORK, double* RWORK, int*   INFO );
    
    // eigenvalue decomposition (symmetric matrix)
    extern "C" void ssyev_( char* JOBZ,  char* UPLO,  int* N, float* Matrix,  int* LDA, float* eigenvalues, float* WORK,  int* LWORK, int* INFO );
    extern "C" void dsyev_( char* JOBZ,  char* UPLO,  int* N, double* Matrix,  int* LDA, double* eigenvalues, double* WORK,  int* LWORK, int* INFO );
    extern "C" void cheev_( char* jobz,  char* uplo,  int* n, std::complex<float>* a,  int* lda, float* w, std::complex<float>* work,  int* lwork, float* rwork, int* info );
    extern "C" void zheev_( char* jobz,  char* uplo,  int* n, std::complex<double>* a,  int* lda, double* w, std::complex<double>* work,  int* lwork, double* rwork, int* info );
    
    // eigenvalue decomposition with conditions (symmetric matrix)
    extern "C" void ssyevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, float * A, int* LDA, 
			     float* VL, float* VU, int* IL, int* UL, float* ABSTOL, int* M,
			     float* W, float* Z, int* LDZ, 
			     float* WORK,  int* LWORK, int* IWORK, int* IFAIL, int* INFO );
    
    extern "C" void dsyevx_( char* JOBZ,  char* RANGE, char* UPLO,  int* N, double* A, int* LDA, 
			     double* VL, double* VU, int* IL, int* UL, double* ABSTOL, int* M,
			     double* W, double* Z, int* LDZ, 
			     double* WORK,  int* LWORK, int* IWORK, int* IFAIL, int* INFO );
    /*
      extern "C" void cheevx_( char* jobz,  char* uplo,  int* n, std::complex<float>* a,  int* lda, 
      float* w, std::complex<float>* work,  int* lwork, float* rwork, int* info );
      
      extern "C" void zheevx_( char* jobz,  char* uplo,  int* n, std::complex<double>* a,  int* lda,
      double* w, std::complex<double>* work,  int* lwork, double* rwork, int* info );
    */
    
    // SVD decomposition
    extern "C" void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A              , int* LDA, double* S, double* U              , int* LDU, double* VT              , int* LDVT, double* WORK, int* LWORK, int* INFO);
    extern "C" void zgesvd_(char* JOBU, char* JOBVT, int* M, int* N, std::complex<double>* A, int* LDA, double* S, std::complex<double>* U, int* LDU, std::complex<double>* VT, int* LDVT, std::complex<double>* WORK, int* LWORK, double* RWORK, int* INFO );
    
    
    // least square fit
    extern "C" void sgelss_(int* M, int* N, int* NRHS,              float*   A, int* LDA,              float*   B, int* LDB, float*  S, float*  RCOND, int* RANK,              float*   WORK, int* LWORK, int* INFO);
    extern "C" void dgelss_(int* M, int* N, int* NRHS,              double*  A, int* LDA,              double*  B, int* LDB, double* S, double* RCOND, int* RANK,              double*  WORK, int* LWORK, int* INFO);
    extern "C" void cgelss_(int* M, int* N, int* NRHS, std::complex<float>*  A, int* LDA, std::complex<float>*  B, int* LDB, float*  S, float*  RCOND, int* RANK, std::complex<float>*  WORK, int* LWORK, float*  RWORK, int* INFO);
    extern "C" void zgelss_(int* M, int* N, int* NRHS, std::complex<double>* A, int* LDA, std::complex<double>* B, int* LDB, double* S, double* RCOND, int* RANK, std::complex<double>* WORK, int* LWORK, double* RWORK, int* INFO);

    // generalized least square fit
    extern "C" void sgglse_(int* M, int* N, int* P, float*  A, int* LDA, float*  B, int* LDB, float*  C, float*  D, float*  X, float* WORK, int* LWORK, int* INFO);
    extern "C" void dgglse_(int* M, int* N, int* P, double* A, int* LDA, double* B, int* LDB, double* C, double* D, double* X, double* WORK, int* LWORK, int* INFO);
    
    // QR decomposition
    extern "C" void sgeqp3_(int* M, int* N, float* A, int* LDA, int* JPVT, float* TAU, float* WORK, int* LWORK, int* INFO);
    extern "C" void dgeqp3_(int* M, int* N, double* A, int* LDA, int* JPVT, double* TAU, double* WORK, int* LWORK, int* INFO);
    
    extern "C" void sorgqr_(int* M, int* N, int* K, float*  A, int* LDA, float* TAU,  float* WORK,  int* LWORK, int* INFO);
    extern "C" void dorgqr_(int* M, int* N, int* K, double* A, int* LDA, double* TAU, double* WORK, int* LWORK, int* INFO);
    
  }
}
    
#endif
