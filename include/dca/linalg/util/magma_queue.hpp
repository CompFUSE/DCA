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

#include <cublas_v2.h>
#include <cuda.h>
#include <cusparse_v2.h>
#include <magma_v2.h>

#include "dca/linalg/util/cuda_stream.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

class MagmaQueue : public linalg::util::CudaStream {
public:
  MagmaQueue() {
    cublasCreate(&cublas_handle_);
    cusparseCreate(&cusparse_handle_);
    int device;
    cudaGetDevice(&device);
    magma_queue_create_from_cuda(device, *this, cublas_handle_, cusparse_handle_, &queue_);
  }

  MagmaQueue(const MagmaQueue& rhs) = delete;

  MagmaQueue(MagmaQueue&& rhs) = default;

  ~MagmaQueue() {
    magma_queue_destroy(queue_);
    cublasDestroy(cublas_handle_);
    cusparseDestroy(cusparse_handle_);
  }

  operator magma_queue_t() const {
    return queue_;
  }

private:
  magma_queue_t queue_ = nullptr;
  cublasHandle_t cublas_handle_ = nullptr;
  cusparseHandle_t cusparse_handle_ = nullptr;
};

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_MAGMA_QUEUE_HPP
