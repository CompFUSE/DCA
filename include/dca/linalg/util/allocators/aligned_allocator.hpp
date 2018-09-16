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

#ifndef DCA_LINALG_UTIL_ALLOCATORS_ALIGNED_ALLOCATOR_HPP
#define DCA_LINALG_UTIL_ALLOCATORS_ALIGNED_ALLOCATOR_HPP

#include <vector>
#include <cuda_runtime.h>

#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename T>
class AlignedAllocator : public std::allocator<T> {
public:
  T* allocate(std::size_t n) {
    T* ptr;
    int err = posix_memalign((void**)&ptr, 128, n * sizeof(T));
    if (err)
      throw(std::bad_alloc());
    return ptr;
  }

  void deallocate(T*& ptr, std::size_t /*n*/ = 0) noexcept {
    free(ptr);
    ptr = nullptr;
  }
};

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_ALIGNED_ALLOCATOR_HPP
