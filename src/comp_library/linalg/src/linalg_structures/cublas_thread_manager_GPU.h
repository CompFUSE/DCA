//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef LIN_ALG_CUBLAS_MANAGER_GPU_H
#define LIN_ALG_CUBLAS_MANAGER_GPU_H

namespace LIN_ALG {

  void create_cublas_handle(int thread_id);
  void destroy_cublas_handle(int thread_id);

  void create_stream_handle(int thread_id);
  void destroy_stream_handle(int thread_id);

 
  void link_thread_stream_to_cublas_handle(int thread_id, int stream_id);
  void synchronize_stream_handle          (int thread_id, int stream_id);

  template<>
  class CUBLAS_THREAD_MANAGER<LIN_ALG::GPU>
  {
  public:

    static void initialize(int id) 
    {
      create_cublas_handle(id);
      create_stream_handle(id);

      link_thread_stream_to_cublas_handle(id, 0);
    }; 

    static void finalize(int id) 
    {
      destroy_cublas_handle(id);
      destroy_stream_handle(id);
    }; 

    static void synchronize_streams(int thread_id, int stream_id)
    {
      synchronize_stream_handle(thread_id, stream_id);
    }

  };
  
}

#endif
