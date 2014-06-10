//-*-C++-*-

#ifndef COPY_FROM_CPU_CPU_H
#define COPY_FROM_CPU_CPU_H

namespace LIN_ALG {

  /*!
   *  \brief structure to copy a matrix from the CPU to the CPU
   */
    template<>
    class COPY_FROM<CPU, CPU>
    {
    public:
      
      template<typename scalartype>
      static void execute(scalartype* source_ptr, scalartype* target_ptr, int size, int thread_id=0, int stream_id=0)
      {
	assert(thread_id>-1 and stream_id>-1);
	
	memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
      }

//       template<typename scalartype>
//       static void execute(scalartype* source_ptr, scalartype* target_ptr, int size, int thread_id, int stream_id){
// 	assert(thread_id>-1 and stream_id>-1);
// 	memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
//       }

//       template<typename scalartype>
//       static void execute(scalartype* source_ptr, scalartype* target_ptr, int size, copy_concurrency_type copy_t=SYNCHRONOUS){
// 	memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
//       }

      template<typename cpu_matrix_type>
      static void execute(cpu_matrix_type& source_cpu_matrix, cpu_matrix_type& target_cpu_matrix, int thread_id=0, int stream_id=0)
      {
	assert(thread_id>-1 and stream_id>-1);

	if(source_cpu_matrix.get_global_size() == target_cpu_matrix.get_global_size())
	  {
	    int size = source_cpu_matrix.get_global_size().first*source_cpu_matrix.get_global_size().second;
	    
	    COPY_FROM<CPU, CPU>::execute(source_cpu_matrix.get_ptr(), target_cpu_matrix.get_ptr(), size);
	  }
	else
	  {
	    COPY_FROM<CPU, CPU>::execute(source_cpu_matrix.get_ptr(), source_cpu_matrix.get_current_size(), source_cpu_matrix.get_global_size(),
					 target_cpu_matrix.get_ptr(), target_cpu_matrix.get_current_size(), target_cpu_matrix.get_global_size());
	  }
      }

      template<typename scalartype>
      static void execute(scalartype* source_ptr, std::pair<int, int>& source_c_s, std::pair<int, int>& source_g_s,
			  scalartype* target_ptr, std::pair<int, int>& target_c_s, std::pair<int, int>& target_g_s,
			  int thread_id=0, int stream_id=0)
      {
	assert(source_c_s == target_c_s);

	assert(thread_id>-1 and stream_id>-1);
	
	for(int j=0; j<source_c_s.second; ++j)
	  memcpy(target_ptr+j*target_g_s.first, source_ptr+j*source_g_s.first, source_c_s.first*sizeof(scalartype));
      }
   };

}

#endif
