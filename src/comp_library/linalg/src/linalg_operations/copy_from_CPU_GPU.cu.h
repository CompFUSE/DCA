//-*-C++-*-

#ifndef COPY_FROM_CPU_TO_GPU_CU_H
#define COPY_FROM_CPU_TO_GPU_CU_H

namespace LIN_ALG 
{
  namespace COPY_FROM_CPU_to_GPU 
  {
    /***************************
     ***   memcopy_h_to_d    ***
     ***************************/

    template<typename scalartype>
    void memcopy_h_to_d(scalartype* target_ptr, scalartype* source_ptr, int size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaMemcpy(target_ptr, source_ptr, sizeof(scalartype)*size, cudaMemcpyHostToDevice);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMemcpyHostToDevice error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();
	  
	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void memcopy_h_to_d(bool*   target_ptr, bool*   source_ptr, int size);
    template void memcopy_h_to_d(int*    target_ptr, int*    source_ptr, int size);
    template void memcopy_h_to_d(float*  target_ptr, float*  source_ptr, int size);
    template void memcopy_h_to_d(double* target_ptr, double* source_ptr, int size);

    /*********************************
     ***   memcopy_h_to_d_async    ***
     *********************************/

    template<typename scalartype>
    void memcopy_h_to_d_async(scalartype* target_ptr, scalartype* source_ptr, int size, int thread_id, int stream_id)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaStream_t& stream_handle = get_stream_handle(thread_id, stream_id);

      cudaError_t ret = cudaMemcpyAsync(target_ptr, source_ptr, sizeof(scalartype)*size, cudaMemcpyHostToDevice, stream_handle);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n\t cudaMemcpyHostToDevice Async error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void memcopy_h_to_d_async(bool*   target_ptr, bool*   source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_h_to_d_async(int*    target_ptr, int*    source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_h_to_d_async(float*  target_ptr, float*  source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_h_to_d_async(double* target_ptr, double* source_ptr, int size, int thread_id, int stream_id);

    /******************************
     ***   memcopy_2D_h_to_d    ***
     ******************************/

    template<typename scalartype>
    void memcopy_2D_h_to_d(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			   scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaMemcpy2D(target_ptr, target_g_s.first*sizeof(scalartype),
				     source_ptr, source_g_s.first*sizeof(scalartype), source_c_s.first*sizeof(scalartype), source_c_s.second,
				     cudaMemcpyHostToDevice);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n\t cudaMemcpy2DHostToDevice error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    template void memcopy_2D_h_to_d(bool*   source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    bool*   target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_h_to_d(int*    source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    int*    target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_h_to_d(float*  source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    float*  target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_h_to_d(double* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    double* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);

    /************************************
     ***   memcopy_2D_h_to_d_async    ***
     ************************************/

    template<typename scalartype>
    void memcopy_2D_h_to_d_async(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
				 scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s,
				 int thread_id, int stream_id)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaStream_t& stream_handle = get_stream_handle(thread_id, stream_id);

      cudaError_t ret = cudaMemcpy2DAsync(target_ptr, target_g_s.first*sizeof(scalartype),
					  source_ptr, source_g_s.first*sizeof(scalartype), source_c_s.first*sizeof(scalartype), source_c_s.second,
					  cudaMemcpyHostToDevice, stream_handle);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n\t cudaMemcpy2DHostToDevice (Async) error :: ret code: " << cudaGetErrorString(ret) << std::endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    template void memcopy_2D_h_to_d_async(bool*   source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
					  bool*   target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_h_to_d_async(int*    source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  int*    target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_h_to_d_async(float*  source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  float*  target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_h_to_d_async(double* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  double* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);

  }
}

#endif
