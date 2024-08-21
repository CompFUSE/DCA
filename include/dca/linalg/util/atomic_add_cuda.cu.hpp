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

#include "dca/config/haves_defines.hpp"
#include "dca/platform/gpu_definitions.h"
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_complex.h"

namespace dca {
namespace linalg {
// dca::linalg::


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

#elif defined(DCA_HAVE_HIP)
// HIP seems to have some horrible problem with concurrent atomic operations.
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

__device__ double inline atomicAddImpl(float* address, const float val) {
  unsigned long int* address_as_int = (unsigned long int*)address;
  unsigned long int old = *address_as_int, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_int, assumed,
                    __float_as_int(val + __int_as_float(assumed)));
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) }
  } while (assumed != old);
  return __int_as_float(old);
}
  
__device__ void inline atomicAdd(float* address, const float val) {
  atomicAddImpl(address, val);
}

__device__ void inline atomicAdd(double* address, const double val) {
  atomicAddImpl(address, val);
}

__device__ void inline atomicAdd(cuDoubleComplex* address, cuDoubleComplex val) {
  double* a_d = reinterpret_cast<double*>(address);
  atomicAddImpl(a_d, val.x);
  atomicAddImpl(a_d + 1, val.y);
}

  __device__ void inline atomicAdd(magmaFloatComplex* const address, magmaFloatComplex val) {
  double* a_d = reinterpret_cast<double*>(address);
  atomicAddImpl(a_d, val.x);
  atomicAddImpl(a_d + 1, val.y);
}

#else
__device__ void inline atomicAdd(double* address, double val) {
  ::atomicAdd(address, val);
}

__device__ void inline atomicAdd(cuComplex* address, cuComplex val) {
  float* a_f = reinterpret_cast<float*>(address);
  ::atomicAdd(a_f, val.x);
  ::atomicAdd(a_f + 1, val.y);
}

__device__ void inline atomicAdd(float* address, float val) {
  ::atomicAdd(address, val);
}

__device__ void inline atomicAdd(cuDoubleComplex* address, cuDoubleComplex val) {
  double* a_d = reinterpret_cast<double*>(address);
  atomicAdd(a_d, val.x);
  atomicAdd(a_d + 1, val.y);
}
#endif  // atomic operation help

}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ATOMIC_ADD_CUDA_CU_HPP
