//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSM_GPU_H
#define LINALG_TRSM_GPU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_TRSM {

    void dtrsm(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB);
    void dtrsm(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB, int id);
  }

  template<>
  class TRSM<GPU>
  {
  public:

    template<typename scalartype>
    inline static void execute(char uplo, char diag, matrix<scalartype, GPU>& A, matrix<scalartype, GPU>& X,
			       int thread_id, int stream_id)
    {      
      assert(uplo=='U' or uplo=='L');
      assert(diag=='U' or diag=='N');
      
      const scalartype ONE(1);
      
      execute('L', uplo, 'N', diag, X.get_number_of_rows(), X.get_number_of_cols(), ONE, A.get_ptr(), A.get_leading_dimension(), X.get_ptr(), X.get_leading_dimension(), thread_id, stream_id);
    }

    /*
    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, float ALPHA, float* A, int LDA, float* B, int LDB){
      GPU_KERNEL_TRSM::strsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
    }
    */
    
    inline static void execute(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB, 
			       int thread_id, int /*stream_id*/)
    {
      // assert(stream_id==0);
      GPU_KERNEL_TRSM::dtrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB, thread_id);
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
