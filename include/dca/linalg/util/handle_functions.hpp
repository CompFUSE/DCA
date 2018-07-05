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
// This file provides cublas related utilities to cublasHandels.

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

// Singleton handle container.
// If not initialized contains the handle for one thread.
inline std::vector<CublasHandle>& getHandleContainer() {
  static std::vector<CublasHandle> handle_container(1);
  return handle_container;
}

// Creates max_threads cublas handles and at least max_threads *
// StreamContainer::streams_per_thread_ cuda streams.
inline void initializeHandleContainer(const int max_threads) {
  if (getStreamContainer().get_max_threads() < max_threads)
    initializeStreamContainer(max_threads);

  auto& handle_container = getHandleContainer();
  handle_container.resize(max_threads);
  // Set default stream.
  for (int id = 0; id < handle_container.size(); ++id)
    handle_container[id].setStream(getStream(id, 0));
}

// Returns the handle associated with thread 'thread_id'.
// Preconditions: 0 <= thread_id < max_threads.
inline cublasHandle_t getHandle(int thread_id) {
  assert(thread_id >= 0 && thread_id < getHandleContainer().size());
  return getHandleContainer()[thread_id];
}

// Sets the stream of the handle associated with thread 'thread_id'
// with the stream returned by getStream(thread_id, stream_id).
// It returns the handle.
// Preconditions: 0 <= thread_id < max_threads,
//                0 <= stream_id < StreamContainer::get_streams_per_thread().
inline cublasHandle_t getHandle(int thread_id, int stream_id) {
  assert(thread_id >= 0 && thread_id < getHandleContainer().size());
  getHandleContainer()[thread_id].setStream(getStream(thread_id, stream_id));
  return getHandleContainer()[thread_id];
}

#else

// Implement SFINAE.
inline void initializeHandleContainer(int /*max_threads*/) {}

#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
