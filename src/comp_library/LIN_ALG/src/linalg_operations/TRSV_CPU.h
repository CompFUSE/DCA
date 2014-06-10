//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSV_CPU_H
#define LINALG_TRSV_CPU_H

namespace LIN_ALG {

  template<>
  class TRSV<CPU>
  {
  public:
    
    template<typename scalartype>
    static void execute(char uplo, matrix<scalartype, CPU>& A, matrix<scalartype, CPU>& X){
      throw std::logic_error(__FUNCTION__);
    }
    
    inline static void execute(char uplo, char trans, char diag, int n, float* A, int LDA, float* X, int incx){
      BLAS::strsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
    }
    
    inline static void execute(char uplo, char trans, char diag, int n, double* A, int LDA, double* X, int incx){
      BLAS::dtrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
    }
    
    inline static void execute(char uplo, char trans, char diag, int n, std::complex<float>* A, int LDA, std::complex<float>* X, int incx){
      BLAS::ctrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
    }
    
    inline static void execute(char uplo, char trans, char diag, int n, std::complex<double>* A, int LDA, std::complex<double>* X, int incx){
      BLAS::ztrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
    }
  };
}

#endif
