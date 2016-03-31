//-*-C++-*-

#ifndef COPY_FROM_CPU_GPU_H
#define COPY_FROM_CPU_GPU_H

namespace LIN_ALG {

  namespace COPY_FROM_CPU_to_GPU {

    template<typename scalartype>
    void memcopy_h_to_d(scalartype* target_ptr, scalartype* source_ptr, int size);

    template<typename scalartype>
    void memcopy_h_to_d_async(scalartype* target_ptr, scalartype* source_ptr, int size, int thread_id, int stream_id);

    template<typename scalartype>
    void memcopy_2D_h_to_d(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
			   scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s);

    template<typename scalartype>
    void memcopy_2D_h_to_d_async(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s, 
				 scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, int thread_id, int stream_id);
  }

  /*!
   *  \brief structure to copy a matrix from the CPU to the GPU
   */
  template<>
  class COPY_FROM<CPU, GPU>
  {
  public:

    template<typename scalartype>
    inline static void execute(scalartype* ptr_cpu, scalartype* ptr_gpu, int size)
    {
      COPY_FROM_CPU_to_GPU::memcopy_h_to_d(ptr_gpu, ptr_cpu, size);
    }

    template<typename scalartype>
    inline static void execute(scalartype* ptr_cpu, scalartype* ptr_gpu, int size, int thread_id, int stream_id)
    {
      COPY_FROM_CPU_to_GPU::memcopy_h_to_d_async(ptr_gpu, ptr_cpu, size, thread_id, stream_id);
    }
    

    template<typename scalartype>
    inline static void execute(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			       scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s)
    {
      assert(source_c_s==target_c_s);

      COPY_FROM_CPU_to_GPU::memcopy_2D_h_to_d(source_ptr, source_c_s, source_g_s, target_ptr, target_c_s, target_g_s);
    }

    template<typename scalartype>
    inline static void execute(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			       scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s, 
			       int thread_id, int stream_id)
    {
      assert(source_c_s==target_c_s);

      COPY_FROM_CPU_to_GPU::memcopy_2D_h_to_d_async(source_ptr, source_c_s, source_g_s, target_ptr, target_c_s, target_g_s, thread_id, stream_id);
    }



    template<typename cpu_matrix_type, typename gpu_matrix_type>
    inline static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix)
    {
      if(cpu_matrix.get_global_size() == gpu_matrix.get_global_size())
	{
	  int size = cpu_matrix.get_global_size().first*cpu_matrix.get_global_size().second;
	  
	  execute(cpu_matrix.get_ptr(), gpu_matrix.get_ptr(), size);
	}
      else
	execute(cpu_matrix.get_ptr(), cpu_matrix.get_current_size(), cpu_matrix.get_global_size(),
		gpu_matrix.get_ptr(), gpu_matrix.get_current_size(), gpu_matrix.get_global_size());
    }

    template<typename cpu_matrix_type, typename gpu_matrix_type>
    inline static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix, int thread_id, int stream_id)
    {
      if(cpu_matrix.get_global_size() == gpu_matrix.get_global_size())
	{
	  int size = cpu_matrix.get_global_size().first*cpu_matrix.get_global_size().second;
	  
	  execute(cpu_matrix.get_ptr(), gpu_matrix.get_ptr(), size, thread_id, stream_id);
	}
      else
	execute(cpu_matrix.get_ptr(), cpu_matrix.get_current_size(), cpu_matrix.get_global_size(),
		gpu_matrix.get_ptr(), gpu_matrix.get_current_size(), gpu_matrix.get_global_size(), 
		thread_id, stream_id);
    }

  };

}

#endif
