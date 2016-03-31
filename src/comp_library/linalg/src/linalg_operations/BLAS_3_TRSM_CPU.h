//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSM_CPU_H
#define LINALG_TRSM_CPU_H

namespace LIN_ALG {

  template<>
  class TRSM<CPU>
  {
  public:
    
    template<typename scalartype>
    inline static void execute(char uplo, char diag, matrix<scalartype, CPU>& A, matrix<scalartype, CPU>& X,  int /*thread_id*/, int /*stream_id*/)
    {      
      assert(uplo=='U' or uplo=='L');
      assert(diag=='U' or diag=='N');
      
      execute('L', uplo, 'N', diag, X.get_number_of_rows(), X.get_number_of_cols(), scalartype(1), A.get_ptr(), A.get_leading_dimension(), X.get_ptr(), X.get_leading_dimension());
    }
    
    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, float ALPHA, float* A, int LDA, float* B, int LDB){
      BLAS::strsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
    }

    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB){
      BLAS::dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
    }

    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, std::complex<float> ALPHA, std::complex<float>* A, int LDA, std::complex<float>* B, int LDB){
      BLAS::ctrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
    }

    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, std::complex<double> ALPHA, std::complex<double>* A, int LDA, std::complex<double>* B, int LDB){
      BLAS::ztrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N, &ALPHA, A, &LDA, B, &LDB);
    }
  };
}

#endif
