// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides a cuda stream container

#ifndef DCA_LINALG_UTIL_STREAM_CONTAINER_HPP
#define DCA_LINALG_UTIL_STREAM_CONTAINER_HPP

#include <array>
#include <cassert>
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <size_t max_threads, size_t streams_per_thread>
class StreamContainer {
public:
  StreamContainer() {
    for (size_t i = 0; i < streams_.size(); ++i) {
      cudaError_t ret = cudaStreamCreate(&(streams_[i]));
      checkRC(ret);
    }
  }
  ~StreamContainer() {
    for (size_t i = 0; i < streams_.size(); ++i) {
      cudaError_t ret = cudaStreamDestroy(streams_[i]);
      checkRC(ret);
    }
  }

  StreamContainer(const StreamContainer<max_threads, streams_per_thread>&) = delete;
  StreamContainer& operator=(const StreamContainer<max_threads, streams_per_thread>&) = delete;

  // Returns the 'stream_id'-th stream associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
  //                0 <= stream_id < DCA_STREAMS_PER_THREADS.
  cudaStream_t operator()(int thread_id, int stream_id) {
    assert(thread_id >= 0 && static_cast<size_t>(thread_id) < max_threads);
    assert(stream_id >= 0 && static_cast<size_t>(stream_id) < streams_per_thread);
    return streams_[stream_id + streams_per_thread * thread_id];
  }

  // Synchronizes the 'stream_id'-th stream associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
  //                0 <= stream_id < DCA_STREAMS_PER_THREADS.
  void sync(int thread_id, int stream_id) {
    cudaError_t ret = cudaStreamSynchronize(operator()(thread_id, stream_id));
    checkRC(ret);
  }

private:
  std::array<cudaStream_t, max_threads * streams_per_thread> streams_;
};
}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_STREAM_CONTAINER_HPP
