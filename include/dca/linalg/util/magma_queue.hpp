// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// RAII wrapper for magma queues.

#ifndef DCA_LINALG_UTIL_MAGMA_QUEUE_HPP
#define DCA_LINALG_UTIL_MAGMA_QUEUE_HPP
#ifdef DCA_HAVE_CUDA

#include <cuda.h>
#include <magma.h>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

class MagmaQueue {
public:
  MagmaQueue() noexcept {
    magma_queue_create(&queue_);
  }

  MagmaQueue(const MagmaQueue&) = delete;
  MagmaQueue& operator=(const MagmaQueue&) = delete;

  MagmaQueue(MagmaQueue&& rhs) noexcept {
    std::swap(queue_, rhs.queue_);
  }

  MagmaQueue& operator=(MagmaQueue&& rhs) noexcept {
    std::swap(queue_, rhs.queue_);
    return *this;
  }

  ~MagmaQueue() {
    magma_queue_destroy(queue_);
  }

  operator magma_queue_t() const noexcept {
    return queue_;
  }

  cudaStream_t getStream() const noexcept {
    return magma_queue_get_cuda_stream(queue_);
  }

private:
  magma_queue_t queue_ = nullptr;
};

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_QUEUE_HPP
