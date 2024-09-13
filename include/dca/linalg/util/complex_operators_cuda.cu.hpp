// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides operators (+,-,*,/) for cuComplex and cuDoubleComplex.
 *  For HIP we just use the magma_operators and magmaComplex types
 */
#ifndef DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
#define DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP

#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/gpu_type_mapping.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

// magma_operators.h defines these operators
#ifdef DCA_HAVE_CUDA

template <typename T>
__device__ __host__ inline void assign(T& a, const T b) {
  a = b;
}

__device__ __host__ inline void assign(cuComplex& a, const std::complex<float>& b) {
  a.x = reinterpret_cast<const float(&)[2]>(b)[0];
  a.y = reinterpret_cast<const float(&)[2]>(b)[1];
}

__device__ __host__ inline void assign(std::complex<float>& a, const cuComplex& b) {
  reinterpret_cast<float(&)[2]>(a)[0] = b.x;
  reinterpret_cast<float(&)[2]>(a)[1] = b.y;
}

__device__ __host__ inline void assign(cuDoubleComplex& a, const std::complex<double>& b) {
  a.x = reinterpret_cast<const double(&)[2]>(b)[0];
  a.y = reinterpret_cast<const double(&)[2]>(b)[1];
}

__device__ __host__ inline void assign(std::complex<double>& a, const cuDoubleComplex& b) {
  reinterpret_cast<double(&)[2]>(a)[0] = b.x;
  reinterpret_cast<double(&)[2]>(a)[1] = b.y;
}

// template <typename T>
// __device__ __host__ inline void assign(dca::util::GPUComplex<T>& a, const std::complex<T>& b) {
//   a.x = reinterpret_cast<const T(&)[2]>(b)[0];
//   a.y = reinterpret_cast<const T(&)[2]>(b)[1];
// }

// template <typename T>
// __device__ __host__ inline void assign(std::complex<T>& a, const dca::util::CUDAComplex<T>& b) {
//   reinterpret_cast<T(&)[2]>(a)[0] = b.x;
//   reinterpret_cast<T(&)[2]>(a)[1] = b.y;
// }

// template <typename T>
// __device__ __host__ inline void assign(dca::util::CUDAComplex<T>& a, const int8_t b) {
//   a.x = static_cast<T>(b);
//   a.y = 0.0;
// }

__device__ __host__ static __inline__ cuComplex operator+(const cuComplex a, const cuComplex b) {
  return cuCaddf(a, b);
}
__device__ __host__ static __inline__ cuComplex operator-(const cuComplex a, const cuComplex b) {
  return cuCsubf(a, b);
}
__device__ __host__ static __inline__ cuComplex operator*(const cuComplex a, const cuComplex b) {
  return cuCmulf(a, b);
}
__device__ __host__ static __inline__ cuComplex operator/(const cuComplex a, const cuComplex b) {
  return cuCdivf(a, b);
}

__device__ __host__ static __inline__ cuComplex operator*=(cuComplex& a, const cuComplex b) {
  a = a * b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator+=(cuComplex& a, const cuComplex b) {
  a = a + b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator+=(cuComplex& a, const float b) {
  a.x = a.x + b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator-=(cuComplex& a, const cuComplex b) {
  a = a - b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator-=(cuComplex& a, const float b) {
  a.x = a.x - b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator-(cuComplex a) {
  a.x = -a.x;
  a.y = -a.y;
  return a;
}
__device__ __host__ static __inline__ cuComplex conj(cuComplex a) {
  a.y = -a.y;
  return a;
}
__device__ __host__ static __inline__ bool operator==(cuComplex a, cuComplex b) {
  return (a.x == b.x) && (a.y == b.y);
}

__device__ __host__ static __inline__ cuComplex operator-(const cuComplex a, const float b) {
  cuComplex tmp = a;
  tmp.x -= b;
  return tmp;
}
__device__ __host__ static __inline__ cuComplex operator-(const float a, const cuComplex b) {
  cuComplex tmp = -b;
  tmp += a;
  return tmp;
}
__device__ __host__ static __inline__ cuComplex operator*(const cuComplex a, const float b) {
  cuComplex tmp = a;
  tmp.x = tmp.x * b;
  tmp.y = tmp.y * b;
  return tmp;
}
__device__ __host__ static __inline__ cuComplex operator*(const float a, const cuComplex b) {
  cuComplex tmp = b;
  tmp.x = tmp.x * a;
  tmp.y = tmp.y * a;
  return tmp;
}
__device__ __host__ static __inline__ cuComplex operator/(const float a, const cuComplex b) {
  cuComplex tmp = b;
  float sq_elem = tmp.x * tmp.x + tmp.y * tmp.y;
  tmp.x = (a * tmp.x) / sq_elem;
  tmp.y = (a * tmp.y) / sq_elem;
  return tmp;
}
__device__ __host__ static __inline__ cuComplex operator/(const cuComplex a, const float b) {
  cuComplex tmp = a;
  tmp.x = (tmp.x / b);
  tmp.y = (tmp.y / b);
  return tmp;
}

__device__ __host__ static __inline__ cuDoubleComplex operator+(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCadd(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCsub(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex a,
                                                                const double b) {
  return {a.x - b, a.y};
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCmul(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCdiv(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex b) {
  a = a * b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator+=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex b) {
  a = a + b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator+=(cuDoubleComplex& a, const double b) {
  a.x = a.x + b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator-=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex b) {
  a = a - b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator-=(cuDoubleComplex& a, const double b) {
  a.x = a.x - b;
  return a;
}

__device__ __host__ static __inline__ cuDoubleComplex operator-(cuDoubleComplex a) {
  a.x = -a.x;
  a.y = -a.y;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex a,
                                                                const double b) {
  return make_cuDoubleComplex(a.x * b, a.y * b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const double a,
                                                                const cuDoubleComplex b) {
  return make_cuDoubleComplex(a * b.x, a * b.y);
}
__device__ __host__ static __inline__ cuDoubleComplex conj(cuDoubleComplex a) {
  a.y = -a.y;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const double a,
                                                                const cuDoubleComplex b) {
  cuDoubleComplex tmp = b;
  double sq_elem = tmp.x * tmp.x + tmp.y * tmp.y;
  tmp.x = (a * tmp.x) / sq_elem;
  tmp.y = (a * tmp.y) / sq_elem;
  return tmp;
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const cuDoubleComplex a,
                                                                const double b) {
  cuDoubleComplex tmp = a;
  tmp.x = (tmp.x / b);
  tmp.y = (tmp.y / b);
  return tmp;
}

__device__ __host__ static __inline__ bool operator==(cuDoubleComplex a, cuDoubleComplex b) {
  return (a.x == b.x) && (a.y == b.y);
}

template <typename Real>
struct Real2CudaComplex;

template <>
struct Real2CudaComplex<double> {
  using type = cuDoubleComplex;
};
template <>
struct Real2CudaComplex<float> {
  using type = cuComplex;
};

template <typename Real>
using GPUComplex = typename Real2CudaComplex<Real>::type;
template <typename Real>
using CUDAComplex = typename Real2CudaComplex<Real>::type;

  
template <class T>
__device__ __host__ static __inline__ CUDAComplex<T>& operator/=(CUDAComplex<T>& a, const T& b) {
  a = {a.x / b, a.y / b};
  return a;
}

template <class T>
__device__ __host__ static __inline__ CUDAComplex<T>& operator*=(CUDAComplex<T>& a, const T& b) {
  a = {a.x * b, a.y * b};
  return a;
}

#endif
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
