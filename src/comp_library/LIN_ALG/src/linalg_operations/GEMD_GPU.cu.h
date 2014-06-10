//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEMD_GPU_CU_H
#define LINALG_GEMD_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_GEMD {

    const static int N_threads = 32;

    __global__ void matrix_times_diagonal_kernel(int m, int n, double* M, int LDM, double* D, double* A, int LDA)
    {
      int I = threadIdx.x + blockIdx.x*blockDim.x;

      int I_min =     (blockIdx.x+0)*blockDim.x;
      int I_max = min((blockIdx.x+1)*blockDim.x, m);

      int J_min =     (blockIdx.y+0)*blockDim.x;
      int J_max = min((blockIdx.y+1)*blockDim.x, n);
	
      if(I>=I_min and I<I_max)
	{
	  for(int j=J_min; j<J_max; ++j)
	    A[I+j*LDA] = M[I+j*LDM]*D[j];
	}
    }

    void execute(int m, int n, double* M, int LDM, double* D, double* A, int LDA,
		 int thread_id, int stream_id)
    {
      if(m>0 and n>0)
	{
#ifdef DEBUG_CUDA
	  cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif
	  
	  int grid_x = get_number_of_blocks(m, N_threads);
	  int grid_y = get_number_of_blocks(n, N_threads);
	  
	  dim3 grid (grid_x, grid_y);
	  dim3 block(N_threads);
	  
	  cudaStream_t& stream_handle = LIN_ALG::get_stream_handle(thread_id, stream_id);

	  matrix_times_diagonal_kernel<<<grid, block, 0, stream_handle>>>(m, n, M, LDM, D, A, LDA);

#ifdef DEBUG_CUDA
	  cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
	}
    }

    /*
    __global__ void matrix_times_diagonal_kernel(int m, int n, double* M, int LDM, double* D, double* A, int LDA){

	int I = threadIdx.x + blockIdx.x*blockDim.x;
	int J = threadIdx.y + blockIdx.y*blockDim.y;
	
	if(I<m && J<n)
	    A[I+J*LDA] = M[I+J*LDM]*D[J];
    }

    void execute(int m, int n, double* M, int LDM, double* D, double* A, int LDA){
      
      const static int Nth = get_number_of_threads();
      
      dim3 bl(m/Nth+1  ,n);
      dim3 th(min(m,Nth),1);
      
      matrix_times_diagonal_kernel<<<bl,th>>>(m, n, M, LDM, D, A, LDA);
    }
    */
  }
}

#endif
