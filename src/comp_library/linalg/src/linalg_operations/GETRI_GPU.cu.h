//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRI_GPU_CU_H
#define LINALG_GETRI_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRI {

      int get_work_space(int n){
	  return magma_get_dgetri_nb(n);
      }

      void sgetri(int N, float* A, int LDA, int* IPIV, float* WORK, int LWORK, int* INFO){
	  magma_sgetri_gpu(N, A, LDA, IPIV, WORK, LWORK, INFO);
      }

      void dgetri(int N, double* A, int LDA, int* IPIV, double* WORK, int LWORK, int* INFO){
	  magma_dgetri_gpu(N, A, LDA, IPIV, WORK, LWORK, INFO);
      }
  }
}

#endif
