//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_LASET_GPU_CU_H
#define LINALG_LASET_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_LASET {

    void set_zero(char UPLO, int M, int N, double* A, int LDA)
    {
      // magmablas_dlaset(UPLO, M, N, A, LDA);
      magmablas_dlaset(magma_uplo_const(UPLO), M, N, 0, 0, A, LDA);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    void set_unity(int M, int N, double* A, int LDA)
    {
      // magmablas_dlaset_identity(M, N, A, LDA);
      magmablas_dlaset(magma_uplo_const('A'), M, N, 0, 1, A, LDA);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    const int BLOCK_SIZE_i = 32;
    const int BLOCK_SIZE_j = 32;

    __global__ void set_zero_kernel(int M, int N, double* A, int LDA)
    {
      assert(blockDim.x == BLOCK_SIZE_i);

      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_i; // index of the index-array

      int i_MIN, i_MAX, J_MIN, J_MAX;

      {// min/max of row-index
	i_MIN = BLOCK_SIZE_i*(blockIdx.x+0);
	i_MAX = BLOCK_SIZE_i*(blockIdx.x+1);
	
	i_MIN = max(i_MIN, 0);
	i_MAX = min(i_MAX, M);
      }

      {// min/max of column-index
	J_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
	J_MAX = BLOCK_SIZE_j*(blockIdx.y+1);
	
	J_MIN = max(J_MIN, 0);
	J_MAX = min(J_MAX, N);
      }

      if(I<i_MAX)
	{
	  for(int J=J_MIN; J<J_MAX; ++J)
	    A[I+J*LDA] = 0.;
	}
    }

    void set_zero(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      if(M>0 and N>0)
	{
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

	  int bl_i = dca::util::ceilDiv(M, BLOCK_SIZE_i);
	  int bl_j = dca::util::ceilDiv(N, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_i, bl_j);
	  
	  cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

	  set_zero_kernel<<<blocks, threads, 0, stream_handle>>>(M, N, A, LDA);

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
	}
    }

    __global__ void set_unity_kernel(int M, int N, double* A, int LDA)
    {
      assert(blockDim.x == BLOCK_SIZE_i);

      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_i; // index of the index-array

      int i_MIN, i_MAX, J_MIN, J_MAX;

      {// min/max of row-index
	i_MIN = BLOCK_SIZE_i*(blockIdx.x+0);
	i_MAX = BLOCK_SIZE_i*(blockIdx.x+1);
	
	i_MIN = max(i_MIN, 0);
	i_MAX = min(i_MAX, M);
      }

      {// min/max of column-index
	J_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
	J_MAX = BLOCK_SIZE_j*(blockIdx.y+1);
	
	J_MIN = max(J_MIN, 0);
	J_MAX = min(J_MAX, N);
      }

      if(I<i_MAX)
	{
	  for(int J=J_MIN; J<J_MAX; ++J)
	    A[I+J*LDA] = 0.;

	  if(J_MIN<=I and I<J_MAX)
	    A[I+I*LDA] = 1.;
	}
    }

    void set_unity(int M, int N, double* A, int LDA, int thread_id, int stream_id)
    {
      if(M>0 and N>0)
	{
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

	  int bl_i = dca::util::ceilDiv(M, BLOCK_SIZE_i);
	  int bl_j = dca::util::ceilDiv(N, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_i, bl_j);
	  
	  cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);
	  
	  set_unity_kernel<<<blocks, threads, 0, stream_handle>>>(M, N, A, LDA);

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
	}
    }

  }

}

#endif
