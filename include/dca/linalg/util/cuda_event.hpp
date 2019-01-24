// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// RAII wrapper for cuda events.

#ifndef DCA_LINALG_UTIL_CUDA_EVENT_HPP
#define DCA_LINALG_UTIL_CUDA_EVENT_HPP

#ifdef DCA_HAVE_CUDA
#include <cuda.h>
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA
class CudaEvent {
public:
  CudaEvent() {
    cudaEventCreate(&event_);
  }

  ~CudaEvent() {
    cudaEventDestroy(event_);
  }

  inline operator cudaEvent_t() {
    return event_;
  }

  void record(cudaStream_t stream) {
    cudaEventRecord(event_, stream);
  }

  void block() const {
    cudaEventSynchronize(event_);
  }

  void block(cudaStream_t stream) const {
    cudaStreamWaitEvent(stream, event_, 0);
  }

  operator bool() const {
    return cudaEventQuery(event_);
  }

private:
  cudaEvent_t event_ = nullptr;
};

// Returns the elapsed time in seconds between two recorded events. Blocks host.
float elapsedTime(cudaEvent_t stop, cudaEvent_t start) {
  cudaEventSynchronize(stop);
  float msec(0);
  cudaEventElapsedTime(&msec, start, stop);
  return 1e-3 * msec;
}

#else

// Define a trivial, non-blocking event in case CUDA is not available.
class CudaEvent {
public:
  template <class T>
  void record(const T& /*stream*/) {}

  void block() const {}

  template <class T>
  void block(const T&) const {}

  operator bool() const {
    return true;
  }
};

#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_CUDA_EVENT_HPP
