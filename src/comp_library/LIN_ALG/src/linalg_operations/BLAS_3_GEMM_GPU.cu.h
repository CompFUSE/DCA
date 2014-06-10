//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEMM_GPU_CU_H
#define LINALG_GEMM_GPU_CU_H

namespace LIN_ALG {

   namespace GPU_KERNEL_GEMM {

     void dgemm(char TRANSA, char TRANSB, int M, int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double  beta, double* C, int LDC)
     {
       cublasStatus_t status = cublasDgemm(get_thread_handle(0), 
					   cublas_operation_type(TRANSA), 
					   cublas_operation_type(TRANSB), 
					   M, N, K, &alpha, A, LDA, B, LDB, &beta, C, LDC);
       
       if(status != CUBLAS_STATUS_SUCCESS)
	 cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);
       
#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     void dgemm(char TRANSA, char TRANSB, int M, int N, int K, double alpha, double* A, int LDA, double* B, int LDB, double  beta, double* C, int LDC, int id)
     {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

       cublasStatus_t status = cublasDgemm(get_thread_handle(id), 
					   cublas_operation_type(TRANSA), 
					   cublas_operation_type(TRANSB), 
					   M, N, K, &alpha, A, LDA, B, LDB, &beta, C, LDC);
       
       if(status != CUBLAS_STATUS_SUCCESS)
	 cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);
       
#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

   }

}

#endif
