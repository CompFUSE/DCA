//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef LIN_ALG_CUBLAS_MANAGER_CPU_H
#define LIN_ALG_CUBLAS_MANAGER_CPU_H

namespace LIN_ALG {

  template<>
  class CUBLAS_THREAD_MANAGER<LIN_ALG::CPU>
  {
  public:

    static void initialize(int id) {}; 
    static void finalize  (int id) {}; 

    static void synchronize_streams(int thread_id, int stream_id) {};
  };
  
}

#endif
