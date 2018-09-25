// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a wrapper to cudaAdd for complex types.

#ifndef DCA_LINALG_UTIL_ATOMIC_ADD_CUDA_CU_HPP
#define DCA_LINALG_UTIL_ATOMIC_ADD_CUDA_CU_HPP

#include <cuComplex.h>

namespace dca {
namespace linalg {
// dca::linalg::

__device__ void inline atomicAdd(cuComplex* address, const cuComplex val) {
  float* a_f = reinterpret_cast<float*>(address);
  ::atomicAdd(a_f, val.x);
  ::atomicAdd(a_f + 1, val.y);
}

__device__ void inline atomicAdd(float* address, const float val) {
  ::atomicAdd(address, val);
}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
// Older devices do not have an hardware atomicAdd for double.
// See
// https://stackoverflow.com/questions/12626096/why-has-atomicadd-not-been-implemented-for-doubles
__device__ double inline atomicAddImpl(double* address, const double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) }
  } while (assumed != old);
  return __longlong_as_double(old);
}

__device__ void inline atomicAdd(double* address, const double val) {
  atomicAddImpl(address, val);
}

#else
__device__ void inline atomicAdd(double* address, const double val) {
  ::atomicAdd(address, val);
}
#endif  // __CUDA_ARCH__

__device__ void inline atomicAdd(cuDoubleComplex* address, const cuDoubleComplex val) {
  double* a_d = reinterpret_cast<double*>(address);
  atomicAdd(a_d, val.x);
  atomicAdd(a_d + 1, val.y);
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ATOMIC_ADD_CUDA_CU_HPP
