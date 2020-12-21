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
// This file provides a cuda stream container

#ifndef DCA_LINALG_UTIL_STREAM_CONTAINER_HPP
#define DCA_LINALG_UTIL_STREAM_CONTAINER_HPP

#include <array>
#include <cassert>
#include <functional>
#include <vector>

#include "dca/linalg/util/cuda_stream.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

class StreamContainer {
public:
  StreamContainer(std::size_t max_threads = 0) : streams_(max_threads) {}

  std::size_t get_max_threads() const {
    return streams_.size();
  }

  std::size_t get_streams_per_thread() const {
    return streams_per_thread_;
  }

  void resize(const int max_threads) {
    streams_.resize(max_threads);
  }

  StreamContainer(const StreamContainer&) = delete;
  StreamContainer& operator=(const StreamContainer&) = delete;

  // Returns the 'stream_id'-th stream associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < get_max_threads(),
  //                0 <= stream_id < streams_per_thread_.
  CudaStream& operator()(int thread_id, int stream_id) {
    assert(thread_id >= 0 && thread_id < get_max_threads());
    assert(stream_id >= 0 && stream_id < streams_per_thread_);
    return streams_[thread_id][stream_id];
  }

  // Synchronizes the 'stream_id'-th stream associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < get_max_threads(),
  //                0 <= stream_id < streams_per_thread_.
  void sync(int thread_id, int stream_id) {
    operator()(thread_id, stream_id).sync();
  }

private:
  constexpr static std::size_t streams_per_thread_ = 2;
  std::vector<std::array<CudaStream, streams_per_thread_>> streams_;
};

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_STREAM_CONTAINER_HPP
