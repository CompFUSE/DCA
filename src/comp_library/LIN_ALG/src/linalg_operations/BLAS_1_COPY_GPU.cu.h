//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_COPY_GPU_CU_H
#define LINALG_COPY_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNEL_COPY {

    void dcopy(int length, double* x, int inc_x, double* y, int inc_y)
    {
      cublasStatus_t status = cublasDcopy(get_thread_handle(0), length, x, inc_x, y, inc_y);
      
      if(status != CUBLAS_STATUS_SUCCESS)
	  cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);

#ifdef DEBUG_CUDA
      cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    void dcopy(int length, double* x, int inc_x, double* y, int inc_y, int id)
    {
      cublasStatus_t status = cublasDcopy(get_thread_handle(id), length, x, inc_x, y, inc_y);
      
      if(status != CUBLAS_STATUS_SUCCESS)
	  cublas_error_msg(status, __FUNCTION__, __FILE__, __LINE__);

#ifdef DEBUG_CUDA
      cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    const static int BLOCK_SIZE_x = 32;
    const static int BLOCK_SIZE_y = 8;

    __global__ void many_row_copies_kernel(int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y)
    {
      int index = threadIdx.x + blockIdx.x*BLOCK_SIZE_x;
      
      int l_MIN = BLOCK_SIZE_y*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_y*(blockIdx.y+1);
      
      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_x);
      
      if(index < N_i)
	{
	  int row_index_y = i_y[index];
	  int row_index_x = i_x[index];

	  for(int l=l_MIN; l<l_MAX; ++l)
	    y[ row_index_y + l*inc_y] = x[ row_index_x + l*inc_x];
	}
    }

    void many_row_copies(int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y, int thread_id, int stream_id)
    {
      if(N_i>0 and N_x>0)
	{
	  int bl_x = get_number_of_blocks(N_i, BLOCK_SIZE_x);
	  int bl_y = get_number_of_blocks(N_x, BLOCK_SIZE_y);
	  
	  dim3 threads(BLOCK_SIZE_x);
	  dim3 blocks (bl_x, bl_y);

	  cudaStream_t& stream_handle = LIN_ALG::get_stream_handle(thread_id, stream_id);
	  
	  many_row_copies_kernel<<<blocks, threads, 0, stream_handle>>>(N_x, N_i, i_x, x, inc_x, i_y, y, inc_y);
	}
    }

    __global__ void many_column_copies_kernel(int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y)
    {      
      int I = threadIdx.x + blockIdx.x*BLOCK_SIZE_x;
	
      int l_MIN = BLOCK_SIZE_y*(blockIdx.y+0);
      int l_MAX = BLOCK_SIZE_y*(blockIdx.y+1);
      
      l_MIN = max(l_MIN, 0);
      l_MAX = min(l_MAX, N_i);
	
      if(I<N_x)
	{
	  for(int l=l_MIN; l<l_MAX; ++l)	
	    y[I + i_y[l]*inc_y] = x[I + i_x[l]*inc_x];
	}
    }

    void many_column_copies(int N_x, int N_i, int* i_x, double* x, int inc_x, int* i_y, double* y, int inc_y, int thread_id, int stream_id)
    {
      if(N_i>0 and N_x>0)
	{
	  int bl_x = get_number_of_blocks(N_x, BLOCK_SIZE_x);
	  int bl_y = get_number_of_blocks(N_i, BLOCK_SIZE_y);
	  
	  dim3 threads(BLOCK_SIZE_x);
	  dim3 blocks (bl_x, bl_y);
	  
	  cudaStream_t& stream_handle = LIN_ALG::get_stream_handle(thread_id, stream_id);
	  
	  many_column_copies_kernel<<<blocks, threads, 0, stream_handle>>>(N_x, N_i, i_x, x, inc_x, i_y, y, inc_y);
	}

#ifdef DEBUG_CUDA
      cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
  }
}

#endif
