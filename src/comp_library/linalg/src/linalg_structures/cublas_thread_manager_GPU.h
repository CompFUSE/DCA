//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef LIN_ALG_CUBLAS_MANAGER_GPU_H
#define LIN_ALG_CUBLAS_MANAGER_GPU_H

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/util/handle_functions.hpp"
#include "dca/linalg/util/stream_functions.hpp"

namespace LIN_ALG {

  template<>
  class CUBLAS_THREAD_MANAGER<LIN_ALG::GPU>
  {
  public:

    static void initialize(int /*id*/) {
      dca::linalg::util::getHandleContainer();
    } 

    static void finalize(int /*id*/) {} 

    static void synchronize_streams(int thread_id, int stream_id)
    {
      dca::linalg::util::syncStream(thread_id, stream_id);
    }

  };
  
}

#endif  // DCA_HAVE_CUDA
#endif
