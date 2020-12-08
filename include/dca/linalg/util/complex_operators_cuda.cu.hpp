// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Baluduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides operators (+,-,*,/) for cuComplex and cuDoubleComplex.

#ifndef DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
#define DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP

#include <cuComplex.h>

namespace dca {
namespace linalg {
// dca::linalg::

__device__ __host__ static __inline__ auto makeComplex(float re, float im) {
  return make_cuComplex(re, im);
}
__device__ __host__ static __inline__ auto makeComplex(double re, double im) {
  return make_cuDoubleComplex(re, im);
}

__device__ __host__ static __inline__ auto makeComplex(const cuComplex& x) {
  return x;
}
__device__ __host__ static __inline__ auto makeComplex(const cuDoubleComplex& x) {
  return x;
}
__device__ __host__ static __inline__ auto makeComplex(const float x) {
  return make_cuComplex(x, 0.f);
}
__device__ __host__ static __inline__ auto makeComplex(const double x) {
  return make_cuDoubleComplex(x, 0);
}

template <class T>
__device__ static __inline__ void assign(T& a, const T& b) {
  a = b;
}

template <class T>
__device__ static __inline__ void assign(CudaComplex<T>& a, const T& b) {
  a = makeComplex(b);
}

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
__device__ __host__ static __inline__ cuComplex operator-(cuComplex a) {
  return make_cuComplex(-a.x, -a.y);
}
__device__ __host__ static __inline__ cuComplex operator-(const cuComplex a, const float b) {
  return make_cuComplex(a.x - b, a.y);
}
__device__ __host__ static __inline__ cuComplex operator-(const float a, const cuComplex b) {
  return make_cuComplex(a - b.x, -b.y);
}
__device__ __host__ static __inline__ cuComplex operator*(const cuComplex a, const float b) {
  return make_cuComplex(a.x * b, a.y * b);
}
__device__ __host__ static __inline__ cuComplex operator*(const float a, const cuComplex b) {
  return make_cuComplex(a * b.x, a * b.y);
}
__device__ __host__ static __inline__ cuComplex operator/(float a, const cuComplex& b) {
  return make_cuComplex(a / b.x, a / b.y);
}
__device__ __host__ static __inline__ cuComplex operator+=(cuComplex& a, const cuComplex b) {
  a = a + b;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator-=(cuComplex& a, const cuComplex b) {
  a = a - b;
  return a;
}
__device__ __host__ static __inline__ cuComplex conj(cuComplex a) {
  a.y = -a.y;
  return a;
}

__device__ __host__ static __inline__ cuDoubleComplex operator+(const cuDoubleComplex& a,
                                                                const cuDoubleComplex& b) {
  return cuCadd(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex& a,
                                                                const cuDoubleComplex& b) {
  return cuCsub(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex& a,
                                                                const cuDoubleComplex& b) {
  return cuCmul(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const cuDoubleComplex& a,
                                                                const cuDoubleComplex& b) {
  return cuCdiv(a, b);
}

__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex& a) {
  return make_cuDoubleComplex(-a.x, -a.y);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex a,
                                                                const double b) {
  return make_cuDoubleComplex(a.x - b, a.y);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const float a,
                                                                const cuDoubleComplex b) {
  return make_cuDoubleComplex(a - b.x, -b.y);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex& a,
                                                                const double b) {
  return make_cuDoubleComplex(a.x * b, a.y * b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const double a,
                                                                const cuDoubleComplex& b) {
  return make_cuDoubleComplex(a * b.x, a * b.y);
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(double a, const cuDoubleComplex& b) {
  return make_cuDoubleComplex(a / b.x, a / b.y);
}
__device__ __host__ static __inline__ cuDoubleComplex operator+=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex& b) {
  a = a + b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex operator-=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex& b) {
  a = a - b;
  return a;
}

template <class T>
__device__ __host__ static __inline__ CudaComplex<T>& operator*=(CudaComplex<T>& a,
                                                                 const CudaComplex<T>& b) {
  a = a * b;
  return a;
}

__device__ __host__ static __inline__ cuComplex& operator*=(cuComplex& a, const cuComplex& b) {
  return a = a * b;
}
__device__ __host__ static __inline__ cuDoubleComplex& operator*=(cuDoubleComplex& a,
                                                                  const cuDoubleComplex& b) {
  return a = a * b;
}

template <class T>
__device__ __host__ static __inline__ CudaComplex<T> operator/(CudaComplex<T>& a, const T& b) {
  return makeComplex(a.x / b, a.y / b);
}

template <class T>
__device__ __host__ static __inline__ CudaComplex<T>& operator*=(CudaComplex<T>& a, const T& b) {
  return a = a * b;
}

template <class T>
__device__ __host__ static __inline__ CudaComplex<T>& operator/=(CudaComplex<T>& a, const T& b) {
  return a = a / b;
}

__device__ __host__ static __inline__ cuDoubleComplex conj(const cuDoubleComplex& a) {
  return make_cuDoubleComplex(a.x, -a.y);
}

}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
