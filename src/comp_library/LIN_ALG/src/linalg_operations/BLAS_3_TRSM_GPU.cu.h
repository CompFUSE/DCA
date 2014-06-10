//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_TRSM_GPU_CU_H
#define LINALG_TRSM_GPU_CU_H

namespace LIN_ALG {
  
  namespace GPU_KERNEL_TRSM {
    
    void dtrsm(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB)
    {
#ifdef DEBUG_CUDA
      cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cublasStatus_t status = cublasDtrsm(get_thread_handle(0), 
					  cublas_side_type(SIDE), 
					  cublas_triangle_type(UPLO), 
					  cublas_operation_type(TRANSA), 
					  cublas_diagonal_type(DIAG), 
					  M, N, &ALPHA, A, LDA, B, LDB);

      if(status != CUBLAS_STATUS_SUCCESS)
	cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);

#ifdef DEBUG_CUDA
      cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    void dtrsm(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, double ALPHA, double* A, int LDA, double* B, int LDB, int id)
    {
#ifdef DEBUG_CUDA
      cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cublasStatus_t status = cublasDtrsm(get_thread_handle(id), 
					  cublas_side_type(SIDE), 
					  cublas_triangle_type(UPLO), 
					  cublas_operation_type(TRANSA), 
					  cublas_diagonal_type(DIAG), 
					  M, N, &ALPHA, A, LDA, B, LDB);

      if(status != CUBLAS_STATUS_SUCCESS)
	cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);
      
#ifdef DEBUG_CUDA
      cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
  }

}

#endif
