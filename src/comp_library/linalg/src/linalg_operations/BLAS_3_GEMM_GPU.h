//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEMM_GPU_H
#define LINALG_GEMM_GPU_H

namespace LIN_ALG {

   namespace GPU_KERNEL_GEMM {

     void dgemm(char TRANSA, char TRANSB, int M, int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double  beta, double* C, int LDC);
     void dgemm(char TRANSA, char TRANSB, int M, int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double  beta, double* C, int LDC, int id);
   }

  template<>
  class GEMM<GPU>
  {
  public:

    template<typename scalartype>
    static void execute(scalartype a, matrix<scalartype, GPU>& A, matrix<scalartype, GPU>& B, 
			scalartype b, matrix<scalartype, GPU>& C,
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

      execute(TRANSA, TRANSB, M, N, K, alpha, A.get_ptr(), LDA, B.get_ptr(), LDB, beta, C.get_ptr(), LDC, thread_id, stream_id);
    }
    
    /*
    inline static void execute(char TRANSA, char TRANSB, int M, int N, int K, 
			       float alpha, 
			       float* A, int LDA,
			       float* B, int LDB,
			       float  beta, 
			       float* C, int LDC)
    {
      GPU_KERNEL_GEMM::sgemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
    */

    inline static void execute(char TRANSA, char TRANSB, int M, int N, int K, 
			       double alpha, 
			       double* A, int LDA,
			       double* B, int LDB,
			       double  beta, 
			       double* C, int LDC,
			       int thread_id, int /*stream_id*/)
    {
      // assert(stream_id==0);
      GPU_KERNEL_GEMM::dgemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC, thread_id);
    }

    /*
    inline static void execute(char TRANSA, char TRANSB, int M, int N, int K, 
			       cuComplex alpha, 
			       cuComplex* A, int LDA,
			       cuComplex* B, int LDB,
			       cuComplex  beta, 
			       cuComplex* C, int LDC)
    {
      cublasCgemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
    
    inline static void execute(char TRANSA, char TRANSB, int M, int N, int K, 
			       cuDoubleComplex alpha, 
			       cuDoubleComplex* A, int LDA,
			       cuDoubleComplex* B, int LDB,
			       cuDoubleComplex  beta, 
			       cuDoubleComplex* C, int LDC)
    {
      cublasZgemm(TRANSA, TRANSB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);
    }
    */
  };

}

#endif
