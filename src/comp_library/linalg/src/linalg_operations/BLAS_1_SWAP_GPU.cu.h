//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_SWAP_GPU_CU_H
#define LINALG_SWAP_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL {
	
    const static int BLOCK_SIZE_i = 32; // rows
    const static int BLOCK_SIZE_j = 8;  // cols

    __global__ void swap_rows_kernel(int A_r, int A_c, double* A_ptr, int A_LD, 
				     int N_i, const int* i_s_ptr, const int* i_t_ptr)
    {
      assert(blockDim.x == BLOCK_SIZE_i);

      __shared__ double T[BLOCK_SIZE_i*BLOCK_SIZE_j];

      int t_id    = threadIdx.x;
      int i_array = threadIdx.x + blockIdx.x*BLOCK_SIZE_i; // index of the index-array

      int i_MIN, i_MAX, J_MIN, J_MAX;

      {// min/max of array
	i_MIN = BLOCK_SIZE_i*(blockIdx.x+0);
	i_MAX = BLOCK_SIZE_i*(blockIdx.x+1);
	
	i_MIN = max(i_MIN, 0  );
	i_MAX = min(i_MAX, N_i);
      }

      {// min/max of column-index
	J_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
	J_MAX = BLOCK_SIZE_j*(blockIdx.y+1);
	
	J_MIN = max(J_MIN, 0  );
	J_MAX = min(J_MAX, A_c);
      }
      
      if(i_array<i_MAX)
	{
	  int i_s = i_s_ptr[i_array];
	  int i_t = i_t_ptr[i_array];

	  assert(i_s>-1 and i_s<A_r);
	  assert(i_t>-1 and i_t<A_r);

	  assert(t_id>-1 and t_id<BLOCK_SIZE_i);

	  for(int j=J_MIN; j<J_MAX; ++j)
	    T[t_id+(j-J_MIN)*BLOCK_SIZE_i] = A_ptr[i_s+j*A_LD];

	  for(int j=J_MIN; j<J_MAX; ++j)
	    A_ptr[i_s+j*A_LD] = A_ptr[i_t+j*A_LD];

	  for(int j=J_MIN; j<J_MAX; ++j)
	    A_ptr[i_t+j*A_LD] = T[t_id+(j-J_MIN)*BLOCK_SIZE_i];
	}
    }

    void swap_many_rows(int M_r, int M_c, double* M_ptr, int M_LD, 
			int N_s, const int* i_s_ptr, const int* i_t_ptr,
			int thread_id, int stream_id)
    {
      if(M_r>0 and M_c>0 and N_s>0)
	{
#ifdef DEBUG_CUDA
	  //CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
	  cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

	  int bl_x = get_number_of_blocks(N_s, BLOCK_SIZE_i);
	  int bl_y = get_number_of_blocks(M_c, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_x, bl_y);
	  
	  cudaStream_t& stream_handle = LIN_ALG::get_stream_handle(thread_id, stream_id);

	  swap_rows_kernel<<<blocks, threads, 0, stream_handle>>>(M_r, M_c, M_ptr, M_LD, 
								  N_s, i_s_ptr, i_t_ptr);
	  
#ifdef DEBUG_CUDA
	  //CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
	  cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
	}
    }		     

    __global__ void swap_cols_kernel(int A_r, int A_c, double* A_ptr, int A_LD, 
				     int N_i, const int* i_s_ptr, const int* i_t_ptr)
    {
      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_i;

      int l_MIN = BLOCK_SIZE_j*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_j*(blockIdx.y+1);

      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_i);

      if(I<A_r)
	{
	  for(int l=l_MIN; l<l_MAX; ++l)
	    {
	      assert(i_s_ptr[l]>-1 and i_s_ptr[l]<A_c);
	      assert(i_t_ptr[l]>-1 and i_t_ptr[l]<A_c);

	      int i_j0 = I+i_s_ptr[l]*A_LD;
	      int i_j1 = I+i_t_ptr[l]*A_LD;
	      
	      double t    = A_ptr[i_j0];
	      A_ptr[i_j0] = A_ptr[i_j1];
	      A_ptr[i_j1] = t;	
	    }
	}
    }

    void swap_many_cols(int M_r, int M_c, double* M_ptr, int M_LD, 
			int N_s, const int* i_s_ptr, const int* i_t_ptr,
			int thread_id, int stream_id)
    {
      if(M_r>0 and M_c>0 and N_s>0)
	{
#ifdef DEBUG_CUDA
	  //CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
	  cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

	  int bl_x = get_number_of_blocks(M_r, BLOCK_SIZE_i);
	  int bl_y = get_number_of_blocks(N_s, BLOCK_SIZE_j);
	  
	  dim3 threads(BLOCK_SIZE_i);
	  dim3 blocks (bl_x, bl_y);

	  cudaStream_t& stream_handle = LIN_ALG::get_stream_handle(thread_id, stream_id);
	  
	  swap_cols_kernel<<<blocks, threads, 0, stream_handle>>>(M_r, M_c, M_ptr, M_LD, 
								  N_s, i_s_ptr, i_t_ptr);
	  
#ifdef DEBUG_CUDA
	  //CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id, stream_id);
	  cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
	}
    }		     
    
  }

}

#endif
