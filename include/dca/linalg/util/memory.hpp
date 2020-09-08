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

  // Do nothing for non arithmetic types.
  template <typename ScalarType>
  static std::enable_if_t<std::is_arithmetic<ScalarType>::value == false, void> setToZero(
      ScalarType /*ptr*/, size_t /*size*/) {}
};

#ifdef DCA_HAVE_CUDA
template <>
struct Memory<GPU> {
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

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_MEMORY_HPP
