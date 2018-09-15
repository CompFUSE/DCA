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

#include <cuda.h>
#include <magma.h>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

class CudaStream {
public:
  CudaStream() {
    cudaStreamCreate(&stream_);
  }

  CudaStream(const CudaStream& other) = delete;

  CudaStream(CudaStream&& other) {
    std::swap(stream_, other.stream_);
  }

  ~CudaStream() {
    if (stream_)
      cudaStreamDestroy(stream_);
  }

  operator cudaStream_t() const {
    return stream_;
  }

private:
  cudaStream_t stream_ = nullptr;
};

}  // util
}  // linalg
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_LINALG_UTIL_CUDA_STREAM_HPP
