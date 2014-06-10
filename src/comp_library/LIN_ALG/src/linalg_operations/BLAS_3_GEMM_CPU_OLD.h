//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEMM_CPU_H
#define LINALG_GEMM_CPU_H

namespace LIN_ALG {

  template<>
  class GEMM<CPU>
  {
  public:

    template<typename scalartype>
    static void execute(scalartype a, matrix<scalartype, CPU>& A, matrix<scalartype, CPU>& B,
			scalartype b, matrix<scalartype, CPU>& C,
			int thread_id, int stream_id)
    {
      assert(A.get_current_size().first  == C.get_current_size().first);
      assert(A.get_current_size().second == B.get_current_size().first);
      assert(B.get_current_size().second == C.get_current_size().second);

      char TRANSA = 'N';
      char TRANSB = 'N';

      scalartype alpha = a;
      scalartype beta  = b;

      int M = A.get_current_size().first;
      int K = A.get_current_size().second;
      int N = B.get_current_size().second;

      int LDA = A.get_global_size().first;
      int LDB = B.get_global_size().first;
      int LDC = C.get_global_size().first;

      execute_gemm(TRANSA, TRANSB, M, N, K, alpha, A.get_ptr(), LDA, B.get_ptr(), LDB, beta, C.get_ptr(), LDC);
    }

    template<typename scalartype>
    inline static void execute(char TRANSA, char TRANSB, int M, int N, int K, 
			       scalartype alpha, 
			       scalartype* A, int LDA,
			       scalartype* B, int LDB,
			       scalartype  beta, 
			       scalartype* C, int LDC,
			       int thread_id=0, int stream_id=0)
    {
      execute_gemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
    
  private:

    inline static void execute_gemm(char TRANSA, char TRANSB, int M, int N, int K, 
				    float alpha, 
				    float* A, int LDA,
				    float* B, int LDB,
				    float  beta, 
				    float* C, int LDC)
    {
      BLAS::sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
    }
    
    inline static void execute_gemm(char TRANSA, char TRANSB, int M, int N, int K, 
				    double alpha, 
				    double* A, int LDA,
				    double* B, int LDB,
				    double  beta, 
				    double* C, int LDC)
    {
      BLAS::dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
    }
    
    inline static void execute_gemm(char TRANSA, char TRANSB, int M, int N, int K, 
				    std::complex<float> alpha, 
				    std::complex<float>* A, int LDA,
				    std::complex<float>* B, int LDB,
				    std::complex<float>  beta, 
				    std::complex<float>* C, int LDC)
    {
      BLAS::cgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
    }
    
    inline static void execute_gemm(char TRANSA, char TRANSB, int M, int N, int K, 
				    std::complex<double> alpha, 
				    std::complex<double>* A, int LDA,
				    std::complex<double>* B, int LDB,
				    std::complex<double>  beta, 
				    std::complex<double>* C, int LDC)
    {
      BLAS::zgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
    }
  };

}

#endif
