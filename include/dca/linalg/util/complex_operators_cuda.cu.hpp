// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides operators (+,-,*,/) for cuComplex and cuDoubleComplex.

#ifndef DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
#define DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP

#include <cuComplex.h>

namespace dca {
namespace linalg {
// dca::linalg::

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
  a.x = -a.x;
  a.y = -a.y;
  return a;
}
__device__ __host__ static __inline__ cuComplex operator*(const cuComplex a, const float b) {
  return make_cuComplex(a.x * b, a.y * b);
}
__device__ __host__ static __inline__ cuComplex operator*(const float a, const cuComplex b) {
  return make_cuComplex(a * b.x, a * b.y);
}
__device__ __host__ static __inline__ cuComplex operator*=(cuComplex& a, const cuComplex b) {
  a = a * b;
  return a;
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

__device__ __host__ static __inline__ cuDoubleComplex operator+(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCadd(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCsub(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCmul(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const cuDoubleComplex a,
                                                                const cuDoubleComplex b) {
  return cuCdiv(a, b);
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
__device__ __host__ static __inline__ cuDoubleComplex operator-=(cuDoubleComplex& a,
                                                                 const cuDoubleComplex b) {
  a = a - b;
  return a;
}
__device__ __host__ static __inline__ cuDoubleComplex conj(cuDoubleComplex a) {
  a.y = -a.y;
  return a;
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
