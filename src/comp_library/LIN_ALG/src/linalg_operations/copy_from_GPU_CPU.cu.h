//-*-C++-*-

#ifndef COPY_FROM_GPU_TO_CPU_CU_H
#define COPY_FROM_GPU_TO_CPU_CU_H

namespace LIN_ALG 
{
  namespace COPY_FROM_GPU_to_CPU 
  {
    /***************************
     ***   memcopy_d_to_h    ***
     ***************************/

    template<typename scalartype>
    void memcopy_d_to_h(scalartype* target_ptr, scalartype* source_ptr, int size)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaMemcpy(target_ptr, source_ptr, sizeof(scalartype)*size, cudaMemcpyDeviceToHost);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMemcpyDeviceToHOST error :: ret code: " << cudaGetErrorString(ret) 
	   << " src-ptr : " << source_ptr 
	   << " tgt-ptr : " << target_ptr << endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void memcopy_d_to_h(bool*   target_ptr, bool*   source_ptr, int size);
    template void memcopy_d_to_h(int*    target_ptr, int*    source_ptr, int size);
    template void memcopy_d_to_h(float*  target_ptr, float*  source_ptr, int size);
    template void memcopy_d_to_h(double* target_ptr, double* source_ptr, int size);

    /*********************************
     ***   memcopy_d_to_h_async    ***
     *********************************/

    template<typename scalartype>
    void memcopy_d_to_h_async(scalartype* target_ptr, scalartype* source_ptr, int size, int thread_id, int stream_id)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaMemcpyAsync(target_ptr, source_ptr, sizeof(scalartype)*size, cudaMemcpyDeviceToHost);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMemcpyDeviceToHost error (Async) : " << cudaGetErrorString(ret) 
		  << " src-ptr : " << source_ptr 
		  << " tgt-ptr : " << target_ptr 
		  << endl;
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

    template void memcopy_d_to_h_async(bool*   target_ptr, bool*   source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_d_to_h_async(int*    target_ptr, int*    source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_d_to_h_async(float*  target_ptr, float*  source_ptr, int size, int thread_id, int stream_id);
    template void memcopy_d_to_h_async(double* target_ptr, double* source_ptr, int size, int thread_id, int stream_id);

    /******************************
     ***   memcopy_2D_d_to_h    ***
     ******************************/

    template<typename scalartype>
    void memcopy_2D_d_to_h(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			   scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

      cudaError_t ret = cudaMemcpy2D(target_ptr, target_g_s.first*sizeof(scalartype),
				     source_ptr, source_g_s.first*sizeof(scalartype), source_c_s.first*sizeof(scalartype), source_c_s.second,
				     cudaMemcpyDeviceToHost);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMemcpy2DDeviceToHost error :: ret code: " << cudaGetErrorString(ret) 
	   << " src-ptr : " << source_ptr << "\t [" << source_c_s.first << ", " << source_c_s.second << "]\t [" << source_g_s.first << ", " << source_g_s.second << "]\n" 
	   << " tgt-ptr : " << target_ptr << "\t [" << target_c_s.first << ", " << target_c_s.second << "]\t [" << target_g_s.first << ", " << target_g_s.second << "]\n";
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    template void memcopy_2D_d_to_h(bool*   source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    bool*   target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_d_to_h(int*    source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    int*    target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_d_to_h(float*  source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    float*  target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
    template void memcopy_2D_d_to_h(double* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				    double* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);

    /************************************
     ***   memcopy_2D_d_to_h_async    ***
     ************************************/

    template<typename scalartype>
    void memcopy_2D_d_to_h_async(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
				 scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id)
    {
#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
      
      cudaStream_t& stream_handle = get_stream_handle(thread_id, stream_id);

      cudaError_t ret = cudaMemcpy2DAsync(target_ptr, target_g_s.first*sizeof(scalartype),
					  source_ptr, source_g_s.first*sizeof(scalartype), source_c_s.first*sizeof(scalartype), source_c_s.second,
					  cudaMemcpyDeviceToHost, stream_handle);

      if(ret != cudaSuccess){
	std::stringstream ss;
	ss << "\n cudaMemcpy2DDeviceToHost (Async) error :: ret code: " << cudaGetErrorString(ret)  
	   << " src-ptr : " << source_ptr << "\t [" << source_c_s.first << ", " << source_c_s.second << "]\t [" << source_g_s.first << ", " << source_g_s.second << "]\n" 
	   << " tgt-ptr : " << target_ptr << "\t [" << target_c_s.first << ", " << target_c_s.second << "]\t [" << target_g_s.first << ", " << target_g_s.second << "]\n";
	std::cout << ss.str();

	throw std::logic_error(__FUNCTION__);
      }

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    template void memcopy_2D_d_to_h_async(bool*   source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  bool*   target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_d_to_h_async(int*    source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  int*    target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_d_to_h_async(float*  source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  float*  target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
    template void memcopy_2D_d_to_h_async(double* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
					  double* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);

  }
}

#endif
