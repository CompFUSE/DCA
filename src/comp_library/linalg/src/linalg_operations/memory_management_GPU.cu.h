//-*-C++-*-

#ifndef LIN_ALG_MEMORY_MANAGEMENT_GPU_CU_H
#define LIN_ALG_MEMORY_MANAGEMENT_GPU_CU_H

namespace LIN_ALG {

  namespace MEMORY_MANAGEMENT_ON_GPU {

    /***********************
     ***    set-to zero
     ***********************/
    
    template<typename scalartype>
    __global__ void set_to_zero_kernel(scalartype* A, int LD, int N) 
    { 
      int I = threadIdx.x + blockIdx.x*blockDim.x;

      if(I<N)
	A[I*LD] = static_cast<scalartype>(0); 
    } 
    
    template<typename scalartype>
    void set_to_zero(scalartype* ptr, int LD, int m){
	
      int th = get_number_of_threads();
      int bl = dca::util::ceilDiv(m, th);
      
      set_to_zero_kernel<<<bl, th>>>(ptr, LD, m);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void set_to_zero(bool  * ptr, int LD, int m);
    template void set_to_zero(int   * ptr, int LD, int m);
    template void set_to_zero(float * ptr, int LD, int m);
    template void set_to_zero(double* ptr, int LD, int m);

    template<typename scalartype>
    void set_to_zero(scalartype* ptr, int m){
	
      int th = get_number_of_threads();
      int bl = dca::util::ceilDiv(m, th);
      
      set_to_zero_kernel<<<bl, th>>>(ptr, 1, m);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void set_to_zero(bool  * ptr, int m);
    template void set_to_zero(int   * ptr, int m);
    template void set_to_zero(float * ptr, int m);
    template void set_to_zero(double* ptr, int m);


    /***********************
     ***    remove-row
     ***********************/

    template<typename scalartype>
    __global__ void remove_first_row_kernel(int m, scalartype* A, int LDA) 
    { 
	int col = blockIdx .x; 
	int row = threadIdx.x; 
	
	extern __shared__ scalartype col_data[];

	col_data[row]    = A[1+row + col*LDA];
	__syncthreads();

	A[row + col*LDA] = col_data[row];
    } 

    template<typename scalartype>
    void remove_first_row(int m, int n, scalartype* A, int LDA)
    { 
      dim3 dim_G(n);
      
      if(m<get_device_properties().maxThreadsPerBlock)
	{
	  dim3 dim_B(m-1);
	  remove_first_row_kernel<<<dim_G, dim_B, get_device_properties().maxThreadsPerBlock>>>(m, A, LDA);
	}
      else
	{
	  int Nm      = m/get_device_properties().maxThreadsPerBlock;
	  int delta_m = m - Nm*m;
	  
	  for(int i=0; i<Nm; ++i){
	    dim3 dim_B(get_device_properties().maxThreadsPerBlock);
		remove_first_row_kernel<<<dim_G, dim_B, get_device_properties().maxThreadsPerBlock>>>(get_device_properties().maxThreadsPerBlock, &A[i*get_device_properties().maxThreadsPerBlock], LDA);
	  }

	  if(delta_m>0){
	    dim3 dim_B(delta_m);
	    remove_first_row_kernel<<<dim_G, dim_B, get_device_properties().maxThreadsPerBlock>>>(delta_m, &A[Nm*get_device_properties().maxThreadsPerBlock], LDA);
	  }
	}

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void remove_first_row(int m, int n, double* A, int LDA);

    /***********************
     ***    remove column
     ***********************/

    template<typename scalartype>
    __global__ void remove_first_col_kernel(int m, int n, scalartype* A, int LDA) 
    { 
	int b_id = blockIdx.x; 

	int t_id = threadIdx.x; 
	int n_th = blockDim.x;

	int i = t_id+b_id*n_th;

	for(int j=0; j<n-1; ++j)
	    A[i + j*LDA] = A[i + (j+1)*LDA];
    }
    
    template<typename scalartype>
    void remove_first_col(int m, int n, scalartype* A, int LDA)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      dim3 dim_G(1);	
      
      if(m<get_device_properties().maxThreadsPerBlock)
	{
	  dim3 dim_B(m);
	  remove_first_col_kernel<<<dim_G, dim_B>>>(m, n, A, LDA);
	}
      else
	{
	  int Nm      = m/get_device_properties().maxThreadsPerBlock;
	  int delta_m = m - Nm*m;
	  
	  {
	    dim3 dim_G(Nm);
	    dim3 dim_B(get_device_properties().maxThreadsPerBlock);
	    remove_first_col_kernel<<<dim_G, dim_B>>>(get_device_properties().maxThreadsPerBlock, n, A, LDA);
	  }
	    
	  if(delta_m>0){
	    dim3 dim_G(1);
	    dim3 dim_B(delta_m);
	    remove_first_col_kernel<<<dim_G, dim_B>>>(delta_m, n, &A[Nm*get_device_properties().maxThreadsPerBlock], LDA);
	  }
	}
      
#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    template void remove_first_col(int m, int n, double* A, int LDA);
  }

}

#endif
