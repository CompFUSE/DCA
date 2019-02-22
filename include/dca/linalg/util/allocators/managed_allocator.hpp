// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides allocators with pinned or mapped memory usable with std::vector.

#ifndef DCA_LINALG_UTIL_ALLOCATORS_MANAGED_ALLOCATOR_HPP
#define DCA_LINALG_UTIL_ALLOCATORS_MANAGED_ALLOCATOR_HPP

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA support."
#endif

#include <cuda_runtime.h>

#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename T>
class ManagedAllocator {
protected:
  T* allocate(std::size_t n) {
    if (n == 0)
      return nullptr;

    cudaError_t ret = cudaMallocManaged((void**)&ptr_, n * sizeof(T));
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__,
                        "\t Managed size requested : " + std::to_string(n * sizeof(T)));
      throw(std::bad_alloc());
    }

    if (stream_)
      cudaStreamAttachMemAsync(stream_, ptr_);

    return ptr_;
  }

  void deallocate(T*& ptr, std::size_t /*n*/ = 0) noexcept {
    cudaError_t ret = cudaFree(ptr);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__);
      std::terminate();
    }
    ptr_ = ptr = nullptr;
  }

public:
  void setStream(cudaStream_t stream) {
    stream_ = stream;
    if (ptr_)
      cudaStreamAttachMemAsync(stream, ptr_);
  }

private:
  T* ptr_ = nullptr;
  cudaStream_t stream_ = nullptr;
};

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_MANAGED_ALLOCATOR_HPP
