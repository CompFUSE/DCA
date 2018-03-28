// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides cuda related utilities to cudaStreams.

#ifndef DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP
#define DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP

#ifdef DCA_HAVE_CUDA
#include <cuda_runtime.h>
#include "dca/linalg/util/def.hpp"
#include "dca/linalg/util/stream_container.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA

inline StreamContainer<DCA_MAX_THREADS, DCA_STREAMS_PER_THREAD>& getStreamContainer() {
  static StreamContainer<DCA_MAX_THREADS, DCA_STREAMS_PER_THREAD> stream_container;
  return stream_container;
}

// Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREADS.
inline cudaStream_t getStream(int thread_id, int stream_id) {
  return getStreamContainer()(thread_id, stream_id);
}

// Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREADS.
inline void syncStream(int thread_id, int stream_id) {
  getStreamContainer().sync(thread_id, stream_id);
}

#else

// SFINAE version of syncStream.
inline void syncStream(int /*thread_id*/, int /*stream_id*/) {}

#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_STREAM_FUNCTIONS_HPP
