// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides memory related utility:
// - allocation, deallocation,
// - setToZero.

#ifndef DCA_LINALG_UTIL_MEMORY_HPP
#define DCA_LINALG_UTIL_MEMORY_HPP

#include <cassert>
#include <complex>
#include <cstring>
#include <stdexcept>

#include "dca/linalg/device_type.hpp"
#include "dca/util/ignore.hpp"

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

#ifdef DCA_HAVE_CUDA
    cudaError_t ret = cudaHostAlloc((void**)&ptr, size * sizeof(ScalarType), cudaHostAllocDefault);
    if (ret != cudaSuccess) {
      checkRCMsg(ret,
                 "\t HOST size requested : " + std::to_string(size) + " * " +
                     std::to_string(sizeof(ScalarType)));
      throw(std::bad_alloc());
    }
#else
    int err = posix_memalign((void**)&ptr, 128, size * sizeof(ScalarType));
    if (err)
      throw(std::bad_alloc());
#endif
  }

  template <typename ScalarType>
  static void deallocate(ScalarType*& ptr) {
#ifdef DCA_HAVE_CUDA
    cudaError_t ret = cudaFreeHost(ptr);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__);
      std::terminate();
    }
#else
    free(ptr);
#endif  // DCA_HAVE_CUDA
    ptr = nullptr;
  }

  // Sets the elements to 0. Only defined for arithmetic types and
  // std::complex of aritmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      ScalarType* ptr, size_t size) {
    std::memset(ptr, 0, size * sizeof(ScalarType));
  }
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      std::complex<ScalarType>* ptr, size_t size) {
    std::memset(static_cast<void*>(ptr), 0, size * sizeof(std::complex<ScalarType>));
  }
};

#ifdef DCA_HAVE_CUDA
template <>
struct Memory<GPU> {
  template <typename ScalarType>
  static void allocate(ScalarType*& ptr, size_t size) {
    assert(ptr == nullptr);

    cudaError_t ret = cudaMalloc((void**)&ptr, size * sizeof(ScalarType));
    if (ret != cudaSuccess) {
      checkRCMsg(ret, "\t DEVICE size requested : " + std::to_string(size));
      throw(std::bad_alloc());
    }
  }

  template <typename ScalarType>
  static void deallocate(ScalarType*& ptr) {
    cudaError_t ret = cudaFree(ptr);
    if (ret != cudaSuccess) {
      printErrorMessage(ret, __FUNCTION__, __FILE__, __LINE__);
      std::terminate();
    }

    ptr = nullptr;
  }

  // Sets the elements to 0. Only defined for arithmetic types and
  // std::complex of aritmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      ScalarType* ptr, size_t size) {
    cudaMemset(ptr, 0, size * sizeof(ScalarType));
  }
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == true, void> setToZero(
      std::complex<ScalarType>* ptr, size_t size) {
    cudaMemset(ptr, 0, size * sizeof(std::complex<ScalarType>));
  }

  // Do nothing for non arithmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == false, void> setToZero(
      ScalarType /*ptr*/, size_t /*size*/) {}
};
#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_MEMORY_HPP
