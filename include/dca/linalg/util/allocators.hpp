// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides allocators with pinned or mapped memory usable with std::vector.

#ifndef DCA_LINALG_UTIL_ALLOCATORS_HPP
#define DCA_LINALG_UTIL_ALLOCATORS_HPP

#include <vector>

#ifdef DCA_HAVE_CUDA
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"
#endif

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

#ifdef DCA_HAVE_CUDA
template <typename T>
class PinnedAllocator : public std::allocator<T> {
public:
  using SizeType = std::size_t;
  using Pointer = T*;
  using ConstPointer = const T*;

  template <typename Tp1>
  struct rebind {
    using other = PinnedAllocator<Tp1>;
  };

  Pointer allocate(SizeType n) {
    Pointer address;
    cudaMallocHost(&address, n * sizeof(T));
    return address;
  }

  void deallocate(Pointer& p, SizeType /*n*/ = 0) {
    cudaFreeHost(p);
    p = nullptr;
  }
};

template <typename T>
class ManagedAllocator : public std::allocator<T> {
public:
  using SizeType = std::size_t;
  using Pointer = T*;
  using ConstPointer = const T*;

  template <typename Tp1>
  struct rebind {
    using other = PinnedAllocator<Tp1>;
  };

  Pointer allocate(SizeType n) {
    Pointer address;
    cudaMallocManaged(&address, n * sizeof(T));
    return address;
  }

  void deallocate(Pointer& p, SizeType /*n*/ = 0) {
    cudaFree(p);
    p = nullptr;
  }
};

template <typename T>
using HostVector = std::vector<T, PinnedAllocator<T>>;
template <typename T>
using ManagedVector = std::vector<T, ManagedAllocator<T>>;

#else
template <typename T>
using HostVector = std::vector<T>;
#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_HPP
