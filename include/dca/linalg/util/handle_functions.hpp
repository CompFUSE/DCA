// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cublas related utilities to cublasHandels.

#ifndef DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
#define DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP

#include <cublas_v2.h>
#include "dca/linalg/util/def.hpp"
#include "dca/linalg/util/handle_container.hpp"
#include "dca/linalg/util/stream_functions.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

inline HandleContainer<DCA_MAX_THREADS>& getHandleContainer() {
  static HandleContainer<DCA_MAX_THREADS> handle_container(getStreamContainer());
  return handle_container;
}

// Returns the handle associated with thread 'thread_id'.
// Preconditions: 0 <= thread_id < DCA_MAX_THREADS.
inline cublasHandle_t getHandle(int thread_id) {
  return getHandleContainer()(thread_id);
}

// Sets the stream of the handle associated with thread 'thread_id'
// with the stream returned by getStream(thread_id, stream_id).
// It returns the handle.
// Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREAD.
inline cublasHandle_t getHandle(int thread_id, int stream_id) {
  getHandleContainer().setStream(thread_id, getStream(thread_id, stream_id));
  return getHandleContainer()(thread_id);
}

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_HANDLE_FUNCTIONS_HPP
