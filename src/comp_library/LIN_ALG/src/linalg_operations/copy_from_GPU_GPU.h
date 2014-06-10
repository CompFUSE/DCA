//-*-C++-*-

#ifndef COPY_FROM_GPU_GPU_H
#define COPY_FROM_GPU_GPU_H

namespace LIN_ALG 
{

  namespace COPY_FROM_GPU_to_GPU 
  {

    template<typename scalartype>
    void memcopy_d_to_d(scalartype* target_ptr, scalartype* source_ptr, int size);

    template<typename scalartype>
    void memcopy_d_to_d_async(scalartype* target_ptr, scalartype* source_ptr, int size, int thread_id, int stream_id);

    template<typename scalartype>
    void memcopy_2D_d_to_d(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			   scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);

    template<typename scalartype>
    void memcopy_2D_d_to_d_async(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
				 scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);
  }

  /*!
   *  \brief structure to copy a matrix from the GPU to the GPU
   */
  template<>
  class COPY_FROM<GPU, GPU>
  {
  public:
    
    template<typename scalartype>
    static void execute(scalartype* source_ptr, scalartype* target_ptr, int size)
    {
      //MEMORY_MANAGEMENT_ON_GPU::memcopy_d_to_d(target_ptr, source_ptr, size);
      COPY_FROM_GPU_to_GPU::memcopy_d_to_d(target_ptr, source_ptr, size);
    }
    
    
    template<typename scalartype>
    static void execute(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s)
    {
      //assert(test_sizes(source_c_s, target_c_s));
      COPY_FROM_GPU_to_GPU::memcopy_2D_d_to_d(source_ptr, source_c_s, source_g_s, target_ptr, target_c_s, target_g_s);
    }
    
    template<typename gpu_matrix_type>
    static void execute(gpu_matrix_type& source_gpu_matrix, gpu_matrix_type& target_gpu_matrix)
    {	
      if(source_gpu_matrix.get_global_size() == target_gpu_matrix.get_global_size())
	{
	  int size = source_gpu_matrix.get_global_size().first*source_gpu_matrix.get_global_size().second;
	  
	  COPY_FROM<GPU, GPU>::execute(source_gpu_matrix.get_ptr(), target_gpu_matrix.get_ptr(), size);
	}
      else
	{
	  COPY_FROM<GPU, GPU>::execute(source_gpu_matrix.get_ptr(), source_gpu_matrix.get_current_size(), source_gpu_matrix.get_global_size(),
				       target_gpu_matrix.get_ptr(), target_gpu_matrix.get_current_size(), target_gpu_matrix.get_global_size());
	}
    }

    /*
    static bool test_sizes(std::pair<int, int>& source_c_s,
			   std::pair<int, int>& target_c_s)
    {
      if(source_c_s==target_c_s)
	return true;
      else
	{
	  cout << "\t(" << source_c_s.first << ", " << source_c_s.second << ")\t<==>\t("<< target_c_s.first << ", " << target_c_s.second << ")" << endl;
	  
	  throw std::logic_error(__FUNCTION__);
	}
      
      return false;
    }
    */

  };
  
}

#endif
