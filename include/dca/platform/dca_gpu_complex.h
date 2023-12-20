// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides working vender or magma complex gpu headers.
 */
#ifndef DCA_GPU_COMPLEX_H
#define DCA_GPU_COMPLEX_H
#include <type_traits>
#include <complex>

#if defined(DCA_HAVE_CUDA)
#include <cuComplex.h>
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"

namespace dca {
namespace linalg {

template <typename T>
__device__ __host__ inline void assign(T& a, const T b) {
  a = b;
}
  
template <typename T>
__device__ __host__ inline void assign(CudaComplex<T>& a, const std::complex<T>& b) {
  a.x = reinterpret_cast<const T(&)[2]>(b)[0];
  a.y = reinterpret_cast<const T(&)[2]>(b)[1];
}

template <typename T>
__device__ __host__ inline void assign(std::complex<T>& a, const CudaComplex<T>& b) {
  reinterpret_cast<T(&)[2]>(a)[0] = b.x;
  reinterpret_cast<T(&)[2]>(a)[1] = b.y;
}

template <typename T>
__device__ __host__ inline void assign(CudaComplex<T>& a, const int8_t b) {
  a.x = static_cast<T>(b);
  a.y = 0.0;
}
  
}
}

#elif defined(DCA_HAVE_HIP)
// hipComplex types are faulty so we use the magma complex types and operators
#include <magma_operators.h>
#include "dca/util/cuda2hip.h"

namespace dca {
namespace linalg {

template <typename T>
__device__ __host__ inline void assign(T& a, const T b) {
  a = b;
}

template <typename T>
struct IsMagmaComplex_t
    : std::disjunction<std::is_same<magmaFloatComplex,T>, std::is_same<magmaDoubleComplex, T>, std::false_type> {};

template <typename T>
using IsMagmaComplex = std::enable_if_t<IsMagmaComplex_t<std::decay_t<T>>::value, bool>;

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
 
}
}
#endif

namespace dca {
namespace linalg {

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

}  // namespace linalg
}  // namespace dca


#endif
