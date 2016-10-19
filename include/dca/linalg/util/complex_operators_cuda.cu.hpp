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

__device__ __host__ static __inline__ cuDoubleComplex operator+(const cuDoubleComplex a, const cuDoubleComplex b) {
  return cuCadd(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator-(const cuDoubleComplex a, const cuDoubleComplex b) {
  return cuCsub(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator*(const cuDoubleComplex a, const cuDoubleComplex b) {
  return cuCmul(a, b);
}
__device__ __host__ static __inline__ cuDoubleComplex operator/(const cuDoubleComplex a, const cuDoubleComplex b) {
  return cuCdiv(a, b);
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_COMPLEX_OPERATORS_CUDA_CU_HPP
