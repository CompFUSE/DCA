// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// RAII wrapper for cuda events.

#ifndef DCA_LINALG_UTIL_GPU_EVENT_HPP
#define DCA_LINALG_UTIL_GPU_EVENT_HPP

#include "dca/config/haves_defines.hpp"
#if defined(DCA_HAVE_GPU)
#include "dca/platform/dca_gpu.h"
#endif  // DCA_HAVE_GPU

#include "dca/linalg/util/gpu_stream.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#if defined(DCA_HAVE_GPU)
class GpuEvent {
public:
  GpuEvent() {
    checkRC(cudaEventCreate(&event_));
  }

  ~GpuEvent() {
    checkRC(cudaEventDestroy(event_));
  }

  inline operator cudaEvent_t() {
    return event_;
  }

  void record(const cudaStream_t stream) {
    checkRC(cudaEventRecord(event_, stream));
  }

  void block() const {
    checkRC(cudaEventSynchronize(event_));
  }

  void block(cudaStream_t stream) const {
    checkRC(cudaStreamWaitEvent(stream, event_, 0));
  }

  operator bool() const {
    cudaError_t return_code = cudaEventQuery(event_);
    if (return_code == cudaSuccess)
      return true;
    else
      return false;
  }

private:
  cudaEvent_t event_ = nullptr;
};

// Returns the elapsed time in seconds between two recorded events. Blocks host.
float elapsedTime(cudaEvent_t stop, cudaEvent_t start) {
  checkRC(cudaEventSynchronize(stop));
  float msec(0);
  checkRC(cudaEventElapsedTime(&msec, start, stop));
  return 1e-3 * msec;
}

#else

// Define a trivial, non-blocking event in case GPU is not available.
class GpuEvent {
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

#endif  // DCA_HAVE_GPU

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_GPU_EVENT_HPP
