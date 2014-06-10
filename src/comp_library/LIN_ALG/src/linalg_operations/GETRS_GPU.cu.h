//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GETRS_GPU_CU_H
#define LINALG_GETRS_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GETRS {

    void sgetrs(char TRANS, int N, int NRHS, float* Matrix_A, int LDA, int* IPIV, float* Matrix_B, int LDB, int* INFO){
      magma_sgetrs_gpu(TRANS, N, NRHS, Matrix_A, LDA, IPIV, Matrix_B, LDB, INFO);
    }
    
    void dgetrs(char TRANS, int N, int NRHS, double* Matrix_A, int LDA, int* IPIV, double* Matrix_B, int LDB, int* INFO){
      magma_dgetrs_gpu(TRANS, N, NRHS, Matrix_A, LDA, IPIV, Matrix_B, LDB, INFO);
    }
  } 
}

#endif
