// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a global container providing access to a CUBLAS handle per thread, and
// utilities related to it.

#ifndef DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
#define DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP

#ifdef DCA_HAVE_CUDA
#include <vector>

#include <cublas_v2.h>

#include "dca/linalg/util/cublas_handle.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA

// Global handle container.
inline std::vector<CublasHandle>& getHandleContainer() {
  static std::vector<CublasHandle> handle_container(1);
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
  getHandleContainer()[thread_id].setStream(getStream(thread_id, stream_id));
  return getHandleContainer()[thread_id];
}

#else

// Implement SFINAE.
inline void resizeHandleContainer(int /*max_threads*/) {}

#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
