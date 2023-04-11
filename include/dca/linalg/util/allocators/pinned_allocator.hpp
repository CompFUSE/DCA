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

#include "dca/config/haves_defines.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#else
#error "This file requires GPU support."
#endif 

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename T>
class PinnedAllocator {
public:
  PinnedAllocator() = default;

  using size_type = std::size_t;
  using pointer = T*;
  using const_pointer = const T*;
  using value_type = T;

  T* allocate(std::size_t n, const void* /*hint*/ = nullptr) {
    if (n == 0)
      return nullptr;
    void* ptr;
    cudaError_t ret = cudaHostAlloc(&ptr, n * sizeof(T), cudaHostAllocDefault);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__,
                        "\t HOST size requested : " + std::to_string(n * sizeof(T)));
      throw(std::bad_alloc());
    }
    return static_cast<T*>(ptr);
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

// These are part of the requirements for a C++ Allocator however these are insufficient
// and they do not appear to be relevant to our codebase yet.
// They are part of what's needed to do a std::move on a PinnedAllocator backed object
// and avoid deallocation and reallocation.
template <class T, class U>
bool operator==(const PinnedAllocator<T>&, const PinnedAllocator<U>&) { return true; }
template <class T, class U>
bool operator!=(const PinnedAllocator<T>&, const PinnedAllocator<U>&) { return false; }

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_PINNED_ALLOCATOR_HPP
