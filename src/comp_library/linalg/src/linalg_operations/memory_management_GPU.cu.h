//-*-C++-*-

#ifndef LIN_ALG_MEMORY_MANAGEMENT_GPU_CU_H
#define LIN_ALG_MEMORY_MANAGEMENT_GPU_CU_H

namespace LIN_ALG {

  namespace MEMORY_MANAGEMENT_ON_GPU {

    /***********************
     ***        get
     ***********************/
    
    template<typename scalartype>
    scalartype get(scalartype* ptr)
    {
      double result;
      
      cudaError_t ret = cudaMemcpy(&result, ptr, sizeof(double), cudaMemcpyDeviceToHost);    

      if(ret != cudaSuccess){
	std::cout << "\n get error :: ret code: " << ret 
		  << " ptr : " << ptr 
		  << std::endl;
	
	throw std::logic_error(__FUNCTION__);
      }
      
      return result;
    }

    template bool   get(bool* ptr  );
    template int    get(int* ptr   );
    template float  get(float* ptr );
    template double get(double* ptr);

    template<typename scalartype>
    scalartype get(scalartype* ptr, int index){
      double result;
      
      cudaError_t ret = cudaMemcpy(&result, ptr+index, sizeof(double), cudaMemcpyDeviceToHost);    

      if(ret != cudaSuccess){
	std::cout << "\n get error :: ret code: " << ret 
		  << " ptr : " << ptr 
		  << std::endl;
	
	throw std::logic_error(__FUNCTION__);
      }

      return result;
    }
    
    template bool   get(bool* ptr  ,  int index);
    template int    get(int* ptr   ,  int index);
    template float  get(float* ptr ,  int index);
    template double get(double* ptr,  int index);

    
    /***********************
     ***        set
     ***********************/

    template<typename scalartype>
    void set(scalartype* ptr, scalartype val){
      cudaError_t ret = cudaMemcpy(ptr, &val, sizeof(scalartype), cudaMemcpyHostToDevice);    

      if(ret != cudaSuccess){
	std::cout << "\n set error :: ret code: " << ret 
		  << " ptr : " << ptr 
		  << std::endl;
	
	throw std::logic_error(__FUNCTION__);
      }
    }

    template void set(bool* ptr,  bool val);
    template void set(int* ptr,  int val);
    template void set(float* ptr,  float val);
    template void set(double* ptr,  double val);


    /***********************
     ***        add
     ***********************/

    template<typename scalartype>
    __global__  void add_kernel(scalartype* ptr, scalartype val){
      ptr[0] += val;
    }
    
    template<typename scalartype>
    void add(scalartype* ptr, scalartype val){
        add_kernel<<<1,1>>>(ptr, val);
    }

    template void add(bool* ptr,  bool val);
    template void add(int* ptr,  int val);
    template void add(float* ptr,  float val);
    template void add(double* ptr,  double val);


    /***********************
     ***    alloc and free
     ***********************/
    
    template<typename scalartype>
    void allocate(scalartype*& ptr, int global_size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      assert(ptr==NULL);

      cudaError_t ret = cudaMalloc((void**) &ptr, global_size*sizeof(scalartype));
      
      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMalloc error :: ret code: " << cudaGetErrorString(ret) 
	   << "\t size requested : " << global_size << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template<typename scalartype>
    void allocate_pinned_host_memory(scalartype*& ptr, int global_size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      assert(ptr==NULL);

      cudaError_t ret = cudaHostAlloc((void**) &ptr, global_size*sizeof(scalartype), cudaHostAllocDefault);
      
      if(ret != cudaSuccess){
	std::stringstream ss;	    
	ss << "\n cudaHostAlloc error :: ret code: " << cudaGetErrorString(ret) << "\n";
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template<typename scalartype>
    void allocate(scalartype*& ptr, std::pair<int,int> global_size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      assert(ptr==NULL);

      cudaError_t ret = cudaMalloc((void**) &ptr, global_size.first*global_size.second*sizeof(scalartype));
      
      if(ret != cudaSuccess){	 
	std::stringstream ss;
	ss << "\n cudaMalloc error :: ret code: " << cudaGetErrorString(ret) 
		  << "\t size requested : " << global_size.first << "\t" << global_size.second 
		  << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template<typename scalartype>
    void allocate_pinned_host_memory(scalartype*& ptr, std::pair<int,int> global_size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      assert(ptr==NULL);

      cudaError_t ret = cudaHostAlloc((void**) &ptr, global_size.first*global_size.second*sizeof(scalartype), cudaHostAllocDefault);
      
      if(ret != cudaSuccess){	 
	std::stringstream ss;   
	ss << "\n cudaHostAlloc error :: ret code: " << cudaGetErrorString(ret) << "\n";
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template<typename scalartype>
    void deallocate(scalartype*& ptr)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaFree(ptr);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaFree error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

      ptr=NULL;

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template<typename scalartype>
    void deallocate_pinned_host_memory(scalartype*& ptr)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaFreeHost(ptr);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaFreeHost error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

      ptr=NULL;

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void   allocate                   (bool*& ptr, int global_size);
    template void   allocate_pinned_host_memory(bool*& ptr, int global_size);
    template void   allocate                   (bool*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(bool*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (bool*& ptr);
    template void deallocate_pinned_host_memory(bool*& ptr);

    template void   allocate                   (int*& ptr, int global_size);
    template void   allocate_pinned_host_memory(int*& ptr, int global_size);
    template void   allocate                   (int*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(int*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (int*& ptr);
    template void deallocate_pinned_host_memory(int*& ptr);

    template void   allocate                   (float*& ptr, int global_size);
    template void   allocate_pinned_host_memory(float*& ptr, int global_size);
    template void   allocate                   (float*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(float*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (float*& ptr);
    template void deallocate_pinned_host_memory(float*& ptr);

    template void   allocate                   (double*& ptr, int global_size);
    template void   allocate_pinned_host_memory(double*& ptr, int global_size);
    template void   allocate                   (double*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(double*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (double*& ptr);
    template void deallocate_pinned_host_memory(double*& ptr);

    template void   allocate                   (std::complex<float>*& ptr, int global_size);
    template void   allocate_pinned_host_memory(std::complex<float>*& ptr, int global_size);
    template void   allocate                   (std::complex<float>*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(std::complex<float>*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (std::complex<float>*& ptr);
    template void deallocate_pinned_host_memory(std::complex<float>*& ptr);

    template void   allocate                   (std::complex<double>*& ptr, int global_size);
    template void   allocate_pinned_host_memory(std::complex<double>*& ptr, int global_size);
    template void   allocate                   (std::complex<double>*& ptr, std::pair<int,int> global_size);
    template void   allocate_pinned_host_memory(std::complex<double>*& ptr, std::pair<int,int> global_size);
    template void deallocate                   (std::complex<double>*& ptr);
    template void deallocate_pinned_host_memory(std::complex<double>*& ptr);


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
