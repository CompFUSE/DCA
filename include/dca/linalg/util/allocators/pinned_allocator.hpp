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

#ifndef DCA_LINALG_UTIL_ALLOCATORS_PINNED_ALLOCATOR_HPP
#define DCA_LINALG_UTIL_ALLOCATORS_PINNED_ALLOCATOR_HPP

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
class PinnedAllocator : public std::allocator<T> {
public:
  T* allocate(std::size_t n) {
    if (n == 0)
      return nullptr;
    T* ptr;
    cudaError_t ret = cudaHostAlloc((void**)&ptr, n * sizeof(T), cudaHostAllocDefault);
    if (ret != cudaSuccess) {
      checkRCMsg(
          ret, "\t HOST size requested : " + std::to_string(n) + " * " + std::to_string(sizeof(T)));
      throw(std::bad_alloc());
    }
    return ptr;
  }

  void deallocate(T*& ptr, std::size_t /*n*/ = 0) noexcept {
    cudaError_t ret = cudaFreeHost(ptr);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__);
      std::terminate();
    }
    ptr = nullptr;
  }
};

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_PINNED_ALLOCATOR_HPP
