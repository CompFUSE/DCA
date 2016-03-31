//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSV_GPU_CU_H
#define LINALG_TRSV_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_TRSV {

    /*
    void strsv(char uplo, char trans, char diag, int n, float* A, int LDA, float* X, int incx){
      cublasStrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }

    void dtrsv(char uplo, char trans, char diag, int n, double* A, int LDA, double* X, int incx){
      cublasDtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }

    void ctrsv(char uplo, char trans, char diag, int n, cuComplex* A, int LDA, cuComplex* X, int incx){
      cublasCtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    
    void ztrsv(char uplo, char trans, char diag, int n, cuDoubleComplex* A, int LDA, cuDoubleComplex* X, int incx){
      cublasZtrsv(uplo, trans, diag, n, A, LDA , X, incx);
    }
    */
  }   
}

#endif
