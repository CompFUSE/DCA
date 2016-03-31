//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSV_GPU_H
#define LINALG_TRSV_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_TRSV {

    void strsv(char uplo, char trans, char diag, int n, float* A, int LDA, float* X, int incx);
    void dtrsv(char uplo, char trans, char diag, int n, double* A, int LDA, double* X, int incx);
  }

  template<>
  class TRSV<GPU>
  {
  public:

    static void execute(char uplo, char trans, char diag, int n, float* A, int LDA, float* X, int incx){
      GPU_KERNEL_TRSV::strsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    
    static void execute(char uplo, char trans, char diag, int n, double* A, int LDA, double* X, int incx){
      GPU_KERNEL_TRSV::dtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    
    /*
    static void execute(char uplo, char trans, char diag, int n, cuComplex* A, int LDA, cuComplex* X, int incx){
      cublasCtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    
    static void execute(char uplo, char trans, char diag, int n, cuDoubleComplex* A, int LDA, cuDoubleComplex* X, int incx){
      cublasZtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    */
  };
    
}

#endif
