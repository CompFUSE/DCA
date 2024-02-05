// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides working vender or magma complex gpu headers.
 *  Order of includes is quite brittle for complex types and operators,
 *  check all platforms when making any changes.
 */
#ifndef DCA_GPU_COMPLEX_H
#define DCA_GPU_COMPLEX_H
#include <type_traits>
#include <complex>

#ifdef DCA_HAVE_CUDA
#include <cuComplex.h>
#endif

#ifdef DCA_HAVE_GPU
#include <magma_types.h>
#endif
#include "dca/util/type_fundamentals.hpp"

#if defined(DCA_HAVE_HIP)
// hipComplex types are faulty so we use the magma complex types and operators
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <magma_operators.h>
#pragma GCC diagnostic pop
#include "dca/util/cuda2hip.h"

namespace dca {
namespace linalg {

template <typename T>
struct IsMagmaComplex_t : std::disjunction<std::is_same<magmaFloatComplex, T>,
                                           std::is_same<magmaDoubleComplex, T>, std::false_type> {};

template <typename T>
using IsMagmaComplex = std::enable_if_t<IsMagmaComplex_t<std::decay_t<T>>::value, bool>;

template <typename T>
__device__ __host__ inline void assign(T& a, const T b) {
  a = b;
}

__device__ __host__ inline void assign(magmaFloatComplex& a, const std::complex<float>& b) {
  a.x = reinterpret_cast<const float(&)[2]>(b)[0];
  a.y = reinterpret_cast<const float(&)[2]>(b)[1];
}

__device__ __host__ inline void assign(std::complex<float>& a, const magmaFloatComplex& b) {
  reinterpret_cast<float(&)[2]>(a)[0] = b.x;
  reinterpret_cast<float(&)[2]>(a)[1] = b.y;
}

__device__ __host__ inline void assign(magmaDoubleComplex& a, const std::complex<double>& b) {
  a.x = reinterpret_cast<const double(&)[2]>(b)[0];
  a.y = reinterpret_cast<const double(&)[2]>(b)[1];
}

__device__ __host__ inline void assign(std::complex<double>& a, const magmaDoubleComplex& b) {
  reinterpret_cast<double(&)[2]>(a)[0] = b.x;
  reinterpret_cast<double(&)[2]>(a)[1] = b.y;
}

}  // namespace linalg
}  // namespace dca
#endif

namespace dca {
namespace linalg {
#ifdef DCA_HAVE_GPU
// The contents of the cast come from en.cppreference.com/w/cpp/numeric/complex
template <typename T>
__device__ __host__ inline void assign(std::complex<T>& a, const T b) {
  a = {b, 0.0};
}

__device__ __host__ inline void assign(double2& a, const int8_t b) {
  a.x = static_cast<double>(b);
  a.y = 0.0;
}

__device__ __host__ inline void assign(float2& a, const int8_t b) {
  a.x = static_cast<float>(b);
  a.y = 0.0;
}
#endif
}  // namespace linalg
}  // namespace dca

#endif
