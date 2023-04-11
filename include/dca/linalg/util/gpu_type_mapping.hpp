// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides better type mapping between host and gpu types
 */

#ifndef DCA_GPU_TYPE_MAPPING_HPP
#define DCA_GPU_TYPE_MAPPING_HPP

#include <complex>
#include <type_traits>
#include <string>
#include <memory>
#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu_complex.h"
#include <magma_v2.h>
#include "dca/util/type_mapping.hpp"

namespace dca {
namespace util {

#ifdef DCA_HAVE_GPU

template <typename T>
using CUDATypeMap = typename std::disjunction<
    OnTypesEqual<T, float, float>, OnTypesEqual<T, double, double>, OnTypesEqual<T, float*, float*>,
    OnTypesEqual<T, double*, double*>, OnTypesEqual<T, const float*, const float*>,
    OnTypesEqual<T, const double*, const double*>, OnTypesEqual<T, float**, float**>,
    OnTypesEqual<T, double**, double**>, OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
    OnTypesEqual<T, std::complex<float>, cuComplex>, OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
    OnTypesEqual<T, std::complex<double>*, cuDoubleComplex*>,
    OnTypesEqual<T, std::complex<float>**, cuComplex**>,
    OnTypesEqual<T, std::complex<double>**, cuDoubleComplex**>,
    OnTypesEqual<T, std::complex<float>*, cuComplex*>,
    OnTypesEqual<T, const std::complex<double>*, const cuDoubleComplex*>,
    OnTypesEqual<T, const std::complex<float>*, const cuComplex*>,
    OnTypesEqual<T, const std::complex<double>&, const cuDoubleComplex&>,
    OnTypesEqual<T, const std::complex<float>&, const cuComplex&>,
    OnTypesEqual<T, const std::complex<float>**, const cuComplex**>,
    OnTypesEqual<T, const std::complex<double>**, const cuDoubleComplex**>,
    OnTypesEqual<T, const std::complex<float>* const*, const cuComplex* const*>,
    OnTypesEqual<T, const std::complex<double>* const*, const cuDoubleComplex* const*>,
    default_type<void>>::type;

template <typename T>
__device__ __host__ CUDATypeMap<T> castGPUType(T var) {
  return reinterpret_cast<CUDATypeMap<T>>(var);
}

template <typename T>
using CUDARealAliasMap =
    typename std::disjunction<OnTypesEqual<T, float2, float>, OnTypesEqual<T, double2, double>,
                              OnTypesEqual<T, cuDoubleComplex, double>,
                              OnTypesEqual<T, cuComplex, float>, default_type<void>>::type;

template <typename T>
CUDARealAliasMap<T> realAliasGPU(T var) {
  return reinterpret_cast<CUDARealAliasMap<T>>(var);
}

template <typename T>
struct IsCudaComplex_t : public std::false_type {};

template <>
struct IsCudaComplex_t<cuComplex> : public std::true_type {};

template <>
struct IsCudaComplex_t<cuDoubleComplex> : public std::true_type {};

template <typename T>
using IsCudaComplex = std::enable_if_t<IsCudaComplex_t<T>::value, bool>;

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
using CudaComplex = typename Real2CudaComplex<Real>::type;

#endif
}  // namespace util
}  // namespace dca

#endif
