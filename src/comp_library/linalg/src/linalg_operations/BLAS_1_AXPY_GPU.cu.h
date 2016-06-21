//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                      \

#ifndef LINALG_AXPY_GPU_CU_H
#define LINALG_AXPY_GPU_CU_H

namespace LIN_ALG {

    namespace GPU_KERNEL_AXPY {

      void daxpy(int length, double f, double* a, int inc_a, double* b, int inc_b)
      {
	cublasStatus_t status = cublasDaxpy(get_thread_handle(0), length, &f, a, inc_a, b, inc_b);

	if(status != CUBLAS_STATUS_SUCCESS)
	  cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);

#ifdef DEBUG_CUDA
      cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
      }

      void daxpy(int length, double f, double* a, int inc_a, double* b, int inc_b, int id)
      {
	cublasStatus_t status = cublasDaxpy(get_thread_handle(id), length, &f, a, inc_a, b, inc_b);

	if(status != CUBLAS_STATUS_SUCCESS)
	  cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);

#ifdef DEBUG_CUDA
      cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
      }
      
    }
}

#endif
