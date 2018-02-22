// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
  MagmaQueue() {
    magma_queue_create(&queue_);
  }

  ~MagmaQueue() {
    magma_queue_destroy(queue_);
  }

  inline operator magma_queue_t() {
    return queue_;
  }

  cudaStream_t getStream() const {
    return magma_queue_get_cuda_stream(queue_);
  }

private:
  magma_queue_t queue_ = nullptr;
};

}  // util
}  // linalg
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_QUEUE_HPP
