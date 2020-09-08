// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// RAII wrapper for cuda stream.

#ifndef DCA_LINALG_UTIL_CUDA_STREAM_HPP
#define DCA_LINALG_UTIL_CUDA_STREAM_HPP

#ifdef DCA_HAVE_CUDA
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA

class CudaStream {
public:
  CudaStream() {
    cudaStreamCreate(&stream_);
  }

  CudaStream(const CudaStream& other) = delete;
  CudaStream& operator=(const CudaStream& other) = delete;

  CudaStream(CudaStream&& other) noexcept {
    swap(other);
  }
  CudaStream& operator=(CudaStream&& other) noexcept {
    swap(other);
    return *this;
  }

  void sync() const {
    checkRC(cudaStreamSynchronize(stream_));
  }

  ~CudaStream() {
    if (stream_)
      cudaStreamDestroy(stream_);
  }

  operator cudaStream_t() const {
    return stream_;
  }

  void swap(CudaStream& other) noexcept {
    std::swap(stream_, other.stream_);
  }

private:
  cudaStream_t stream_ = nullptr;
};

#else  // DCA_HAVE_CUDA

// Mock object.
class CudaStream {
public:
  CudaStream() = default;

  void sync() const {}
};

#endif  // DCA_HAVE_CUDA

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_CUDA_STREAM_HPP
