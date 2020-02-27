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
// This file provides access to a global instance of StreamContainer.

#ifndef DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP
#define DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP

#include "dca/linalg/util/stream_container.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

// Global stream container.
inline StreamContainer& getStreamContainer() {
  // Initialize resources for one thread.
  static StreamContainer stream_container(1);
  return stream_container;
}

inline void resizeStreamContainer(const int max_threads) {
  getStreamContainer().resize(max_threads);
}

// Preconditions: 0 <= thread_id < max_threads,
//                0 <= stream_id < StreamContainer::streams_per_thread_.
inline CudaStream& getStream(int thread_id, int stream_id) {
  return getStreamContainer()(thread_id, stream_id);
}

// Preconditions: 0 <= thread_id < max_threads,
//                0 <= stream_id < StreamContainer::streams_per_thread_.
inline void syncStream(int thread_id, int stream_id) {
  getStreamContainer().sync(thread_id, stream_id);
}

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP
