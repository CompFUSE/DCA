// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file provides a global container providing access to a CUBLAS handle per thread, and
// utilities related to it.

#ifndef DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
#define DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP

#include <vector>
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_blas.h"

#include "dca/linalg/util/stream_functions.hpp"
#include "dca/linalg/util/gpuBLAS_handles.hpp"
#include "dca/linalg/util/gpu_stream.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_GPU


// Global handle container.
inline std::vector<GpuBLASHandle>& getHandleContainer() {
  static std::vector<GpuBLASHandle> handle_container(1);
  return handle_container;
}

// Creates max_threads cublas handles and at least max_threads *
// StreamContainer::streams_per_thread_ cuda streams.
inline void resizeHandleContainer(const std::size_t max_threads) {
  if (getStreamContainer().get_max_threads() < max_threads)
    resizeStreamContainer(max_threads);

  getHandleContainer().resize(max_threads);
}

// Returns the handle associated with thread 'thread_id'.
// Preconditions: 0 <= thread_id < max_threads.
inline cublasHandle_t getHandle(const int thread_id) {
  assert(thread_id >= 0 && thread_id < getHandleContainer().size());
  return getHandleContainer()[thread_id];
}

// Returns the handle associated with thread 'thread_id' after setting its cuda stream to the one
// returned by getStream(thread_id, stream_id).
// It returns the handle.
// Preconditions: 0 <= thread_id < max_threads,
//                0 <= stream_id < StreamContainer::get_streams_per_thread().
inline cublasHandle_t getHandle(const int thread_id, const int stream_id) {
  assert(thread_id >= 0 && thread_id < getHandleContainer().size());
  GpuStream& stream = getStream(thread_id, stream_id);
  getHandleContainer()[thread_id].setStream(stream.streamActually());
  return getHandleContainer()[thread_id];
}

#else


inline void resizeHandleContainer(const std::size_t max_threads) {
  if (getStreamContainer().get_max_threads() < max_threads)
    resizeStreamContainer(max_threads);
}

#endif  // DCA_HAVE_GPU

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
