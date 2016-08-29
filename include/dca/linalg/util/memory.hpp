// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file providesmemory related utility:
// - allocation, deallocation.

#ifndef DCA_LINALG_UTIL_MEMORY_HPP
#define DCA_LINALG_UTIL_MEMORY_HPP

#include <cassert>
#include "dca/linalg/device_type.hpp"

#ifdef DCA_HAVE_CUDA
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"
#endif

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <DeviceType device_name>
struct Memory {};

template <>
struct Memory<CPU> {
  template <typename ScalarType>
  static void allocate(ScalarType*& ptr, size_t size) {
    assert(ptr == nullptr);

#if defined(ENABLE_PINNED_MEMORY_ALLOCATION) and defined(DCA_HAVE_CUDA)
    cudaError_t ret = cudaHostAlloc((void**)&ptr, size * sizeof(ScalarType));
    checkRCMsg(ret, "\t size requested : " + std::to_string(size) + " * " +
                        std::to_string(sizeof(ScalarType)));
    checkErrorsCudaDebug();
#else
    posix_memalign((void**)&ptr, 128, size * sizeof(ScalarType));
#endif
  }

  template <typename ScalarType>
  static void deallocate(ScalarType*& ptr) {
#if defined(ENABLE_PINNED_MEMORY_ALLOCATION) and defined(DCA_HAVE_CUDA)
    cudaError_t ret = cudaHostFree(ptr);
    checkRC(ret);
    checkErrorsCudaDebug();
#else
    free(ptr);
#endif
    ptr = nullptr;
  }
};

#ifdef DCA_HAVE_CUDA
template <>
struct Memory<GPU> {
  template <typename ScalarType>
  static void allocate(ScalarType*& ptr, size_t size) {
    assert(ptr == nullptr);

    cudaError_t ret = cudaMalloc((void**)&ptr, size * sizeof(ScalarType));
    checkRCMsg(ret, "\t size requested : " + std::to_string(size));
    checkErrorsCudaDebug();
  }

  template <typename ScalarType>
  static void deallocate(ScalarType*& ptr) {
    cudaError_t ret = cudaFree(ptr);
    checkRC(ret);

    ptr = nullptr;
  }
};
#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_MEMORY_HPP
