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

#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/util/type_mapping.hpp"
#include "dca/util/type_utils.hpp"

#include <magma_v2.h>

namespace dca {
namespace util {

#if defined(DCA_HAVE_CUDA)

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

namespace details {
// dca::linalg::util::details::
template <typename Real>
struct ComplexContainer;
template <>
struct ComplexContainer<double> {
  using type = cuDoubleComplex;
};
template <>
struct ComplexContainer<float> {
  using type = cuComplex;
};
}  // namespace details
// dca::linalg::util::

template <typename Real>
using CudaComplex = typename details::ComplexContainer<Real>::type;

template <typename T>
struct IsCudaComplex_t : public std::false_type {};

template <>
struct IsCudaComplex_t<cuComplex> : public std::true_type {};

template <>
struct IsCudaComplex_t<cuDoubleComplex> : public std::true_type {};

template <typename T, typename = bool>
struct CudaScalar_impl {};

template <typename T>
struct CudaScalar_impl<T, IsReal<T>> {
  using value_type = T;
};

template <typename T>
struct CudaScalar_impl<T, IsComplex<T>> {
  using value_type = CudaComplex<RealAlias<T>>;
};

template <typename T>
using CudaScalar = typename CudaScalar_impl<T>::value_type;

template <typename T>
inline auto GPUTypeConversion(T var, typename std::enable_if_t<IsCudaComplex_t<T>::value>* = 0) {
  return CUDATypeMap<T>{var.x, var.y};
}

template <typename T>
inline auto GPUTypeConversion(
    T var, typename std::enable_if_t<IsComplex_t<T>::value && (!IsCudaComplex_t<T>::value)>* = 0) {
  return CUDATypeMap<T>{var.real(), var.imag()};
}

template <typename T>
inline auto GPUTypeConversion(T var, typename std::enable_if_t<std::is_floating_point<T>::value>* = 0) {
  return var;
}

#endif
}  // namespace util
}  // namespace dca

#endif
