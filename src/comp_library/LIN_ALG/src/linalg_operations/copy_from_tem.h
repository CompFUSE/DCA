//-*-C++-*-

#ifndef COPY_FROM_H
#define COPY_FROM_H

namespace LIN_ALG {

    template<device_type source_device_name, device_type target_device_name>
    class COPY_FROM
    {
      template<typename cpu_matrix_type, typename gpu_matrix_type>
      static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix);

      template<typename cpu_matrix_type, typename gpu_matrix_type>
      static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix, int thread_id, int stream_id);

      template<typename cpu_scalartype, typename gpu_scalartype>
      static void execute(gpu_scalartype*   ptr, std::pair<int, int>&   c_s, std::pair<int, int>&   g_s,
			  cpu_scalartype* o_ptr, std::pair<int, int>& o_c_s, std::pair<int, int>& o_g_s);
    };

}

#endif
