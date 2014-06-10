//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRF_GPU_CU_H
#define LINALG_GETRF_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRF {

      void sgetrf(int M, int N, float* A, int LDA, int* IPIV, int* INFO){
	  magma_sgetrf_gpu(M, N, A, LDA, IPIV, INFO);
      }
      
      void dgetrf(int M, int N, double* A, int LDA, int* IPIV, int* INFO){
	  magma_dgetrf_gpu(M, N, A, LDA, IPIV, INFO);
      }
  }
}

#endif
