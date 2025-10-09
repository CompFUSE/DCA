// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file provides the castCudaComplex utility and the CudaComplex typedef.

#ifndef DCA_LINALG_UTIL_CAST_CUDA_HPP
#define DCA_LINALG_UTIL_CAST_CUDA_HPP

#include <complex>

#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/util/type_mapping.hpp"

#include <magma_v2.h>
namespace dca::linalg::util {

#if defined(DCA_HAVE_CUDA)
// returns a cuComplex pointer.
inline cuComplex** castCudaComplex(std::complex<float>** ptr) {
  return reinterpret_cast<cuComplex**>(ptr);
}
inline cuComplex* castCudaComplex(std::complex<float>* ptr) {
  return reinterpret_cast<cuComplex*>(ptr);
}
inline cuComplex* castCudaComplex(std::complex<float>& el) {
  return castCudaComplex(&el);
}
inline const cuComplex* const* castCudaComplex(const std::complex<float>* const* ptr) {
  return reinterpret_cast<const cuComplex* const*>(ptr);
}
inline const cuComplex* castCudaComplex(const std::complex<float>* ptr) {
  return reinterpret_cast<const cuComplex*>(ptr);
}
inline const cuComplex* castCudaComplex(const std::complex<float>& el) {
  return castCudaComplex(&el);
}

// returns a cuDoubleComplex pointer.
inline cuDoubleComplex** castCudaComplex(std::complex<double>** ptr) {
  return reinterpret_cast<cuDoubleComplex**>(ptr);
}
inline cuDoubleComplex* castCudaComplex(std::complex<double>* ptr) {
  return reinterpret_cast<cuDoubleComplex*>(ptr);
}
inline cuDoubleComplex* castCudaComplex(std::complex<double>& el) {
  return castCudaComplex(&el);
}
inline const cuDoubleComplex* const* castCudaComplex(const std::complex<double>* const* ptr) {
  return reinterpret_cast<const cuDoubleComplex* const*>(ptr);
}
inline const cuDoubleComplex* castCudaComplex(const std::complex<double>* ptr) {
  return reinterpret_cast<const cuDoubleComplex*>(ptr);
}
inline const cuDoubleComplex* castCudaComplex(const std::complex<double>& el) {
  return castCudaComplex(&el);
}
#elif defined(DCA_HAVE_HIP)
// returns a magmaFloatComplex pointer.
inline magmaFloatComplex** castCudaComplex(std::complex<float>** ptr) {
  return reinterpret_cast<magmaFloatComplex**>(ptr);
}
inline magmaFloatComplex* castCudaComplex(std::complex<float>* ptr) {
  return reinterpret_cast<magmaFloatComplex*>(ptr);
}
inline magmaFloatComplex* castCudaComplex(std::complex<float>& el) {
  return castCudaComplex(&el);
}
inline const magmaFloatComplex* const* castCudaComplex(const std::complex<float>* const* ptr) {
  return reinterpret_cast<const magmaFloatComplex* const*>(ptr);
}
inline const magmaFloatComplex* castCudaComplex(const std::complex<float>* ptr) {
  return reinterpret_cast<const magmaFloatComplex*>(ptr);
}
inline const magmaFloatComplex* castCudaComplex(const std::complex<float>& el) {
  return castCudaComplex(&el);
}

// returns a hipDoubleComplex pointer.
inline magmaDoubleComplex** castCudaComplex(std::complex<double>** ptr) {
  return reinterpret_cast<magmaDoubleComplex**>(ptr);
}
inline magmaDoubleComplex* castCudaComplex(std::complex<double>* ptr) {
  return reinterpret_cast<magmaDoubleComplex*>(ptr);
}
inline magmaDoubleComplex* castCudaComplex(std::complex<double>& el) {
  return castCudaComplex(&el);
}
inline const magmaDoubleComplex* const* castCudaComplex(const std::complex<double>* const* ptr) {
  return reinterpret_cast<const magmaDoubleComplex* const*>(ptr);
}
inline const magmaDoubleComplex* castCudaComplex(const std::complex<double>* ptr) {
  return reinterpret_cast<const magmaDoubleComplex*>(ptr);
}
inline const magmaDoubleComplex* castCudaComplex(const std::complex<double>& el) {
  return castCudaComplex(&el);
}
#endif

inline magmaDoubleComplex** castMAGMAComplex(std::complex<double>** ptr) {
  return reinterpret_cast<magmaDoubleComplex**>(ptr);
}
inline magmaDoubleComplex* castMAGMAComplex(std::complex<double>* ptr) {
  return reinterpret_cast<magmaDoubleComplex*>(ptr);
}
inline magmaDoubleComplex* castMAGMAComplex(std::complex<double>& el) {
  return castMAGMAComplex(&el);
}
inline const magmaDoubleComplex* const* castMAGMAComplex(const std::complex<double>* const* ptr) {
  return reinterpret_cast<const magmaDoubleComplex* const*>(ptr);
}
inline const magmaDoubleComplex* castMAGMAComplex(const std::complex<double>* ptr) {
  return reinterpret_cast<const magmaDoubleComplex*>(ptr);
}
inline const magmaDoubleComplex* castMAGMAComplex(const std::complex<double>& el) {
  return castMAGMAComplex(&el);
}
inline magmaFloatComplex** castMAGMAComplex(std::complex<float>** ptr) {
  return reinterpret_cast<magmaFloatComplex**>(ptr);
}
inline magmaFloatComplex* castMAGMAComplex(std::complex<float>* ptr) {
  return reinterpret_cast<magmaFloatComplex*>(ptr);
}
inline magmaFloatComplex* castMAGMAComplex(std::complex<float>& el) {
  return castMAGMAComplex(&el);
}
inline const magmaFloatComplex* const* castMAGMAComplex(const std::complex<float>* const* ptr) {
  return reinterpret_cast<const magmaFloatComplex* const*>(ptr);
}
inline const magmaFloatComplex* castMAGMAComplex(const std::complex<float>* ptr) {
  return reinterpret_cast<const magmaFloatComplex*>(ptr);
}
inline const magmaFloatComplex* castMAGMAComplex(const std::complex<float>& el) {
  return castMAGMAComplex(&el);
}

#ifdef DCA_HAVE_CUDA
#define cublasDoubleComplex cuDoubleComplex
#define cublasComplex cuComplex
#endif

inline cublasDoubleComplex** castCUBLASComplex(std::complex<double>** ptr) {
  return reinterpret_cast<cublasDoubleComplex**>(ptr);
}
inline cublasDoubleComplex* castCUBLASComplex(std::complex<double>* ptr) {
  return reinterpret_cast<cublasDoubleComplex*>(ptr);
}
inline cublasDoubleComplex* castCUBLASComplex(std::complex<double>& el) {
  return castCUBLASComplex(&el);
}
inline const cublasDoubleComplex* const* castCUBLASComplex(const std::complex<double>* const* ptr) {
  return reinterpret_cast<const cublasDoubleComplex* const*>(ptr);
}
inline const cublasDoubleComplex* castCUBLASComplex(const std::complex<double>* ptr) {
  return reinterpret_cast<const cublasDoubleComplex*>(ptr);
}
inline const cublasDoubleComplex* castCUBLASComplex(const std::complex<double>& el) {
  return castCUBLASComplex(&el);
}
inline cublasComplex** castCUBLASComplex(std::complex<float>** ptr) {
  return reinterpret_cast<cublasComplex**>(ptr);
}
inline cublasComplex* castCUBLASComplex(std::complex<float>* ptr) {
  return reinterpret_cast<cublasComplex*>(ptr);
}
inline cublasComplex* castCUBLASComplex(std::complex<float>& el) {
  return castCUBLASComplex(&el);
}
inline const cublasComplex* const* castCUBLASComplex(const std::complex<float>* const* ptr) {
  return reinterpret_cast<const cublasComplex* const*>(ptr);
}
inline const cublasComplex* castCUBLASComplex(const std::complex<float>* ptr) {
  return reinterpret_cast<const cublasComplex*>(ptr);
}
inline const cublasComplex* castCUBLASComplex(const std::complex<float>& el) {
  return castCUBLASComplex(&el);
}

// Provides a templated typedef.
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
template <typename Real>
using CudaComplex = typename details::ComplexContainer<Real>::type;
}  // namespace dca::linalg::util

inline magmaDoubleComplex convertToMagmaType(std::complex<double> var) {
  return {reinterpret_cast<double (&)[2]>(var)[0], reinterpret_cast<double (&)[2]>(var)[1]};
}

inline magmaFloatComplex convertToMagmaType(std::complex<float> var) {
  return {reinterpret_cast<float (&)[2]>(var)[0], reinterpret_cast<float (&)[2]>(var)[1]};
}

#ifdef DCA_HAVE_HIP
inline magmaFloatComplex convertToMagmaType(HIP_vector_type<float, 2> var) {
  return {reinterpret_cast<float (&)[2]>(var)[0], reinterpret_cast<float (&)[2]>(var)[1]};
}

inline magmaDoubleComplex convertToMagmaType(HIP_vector_type<double, 2> var) {
  return {reinterpret_cast<double (&)[2]>(var)[0], reinterpret_cast<double (&)[2]>(var)[1]};
}
#endif

namespace dca::util {
template <typename T>
using MAGMATypeMap = typename std::disjunction<
    OnTypesEqual<T, float, float>, OnTypesEqual<T, double, double>, OnTypesEqual<T, float*, float*>,
    OnTypesEqual<T, double*, double*>, OnTypesEqual<T, const float*, const float*>,
    OnTypesEqual<T, const double*, const double*>, OnTypesEqual<T, float**, float**>,
    OnTypesEqual<T, const float**, const float**>, OnTypesEqual<T, double**, double**>,
    OnTypesEqual<T, const double**, const double**>,
    OnTypesEqual<T, std::complex<double>*, magmaDoubleComplex*>,
    OnTypesEqual<T, std::complex<float>**, magmaFloatComplex**>,
    OnTypesEqual<T, std::complex<double>**, magmaDoubleComplex**>,
    OnTypesEqual<T, std::complex<float>*, magmaFloatComplex*>,
    OnTypesEqual<T, float2, magmaFloatComplex>, OnTypesEqual<T, double2, magmaDoubleComplex>,
    OnTypesEqual<T, const std::complex<double>*, const magmaDoubleComplex*>,
    OnTypesEqual<T, const std::complex<float>*, const magmaFloatComplex*>,
    OnTypesEqual<T, const std::complex<double>&, const magmaDoubleComplex&>,
    OnTypesEqual<T, const std::complex<float>&, const magmaFloatComplex&>,
    OnTypesEqual<T, const std::complex<float>**, const magmaFloatComplex**>,
    OnTypesEqual<T, const std::complex<double>**, const magmaDoubleComplex**>,
    OnTypesEqual<T, const std::complex<float>* const*, const magmaFloatComplex* const*>,
    OnTypesEqual<T, const std::complex<double>* const*, const magmaDoubleComplex* const*>,
    OnTypesEqual<T, const double2* const*, const magmaDoubleComplex* const*>,
    OnTypesEqual<T, const float2* const*, const magmaFloatComplex* const*>,
#ifdef DCA_HAVE_HIP
    OnTypesEqual<T, const HIP_vector_type<float, 2>* const*, const magmaFloatComplex* const*>,
    OnTypesEqual<T, const HIP_vector_type<double, 2>* const*, const magmaDoubleComplex* const*>,
#endif
    default_type<void>>::type;
template <typename T>
__device__ __host__ MAGMATypeMap<T> castMagmaType(T var) {
  return reinterpret_cast<MAGMATypeMap<T>>(var);
}

}  // namespace dca::util

#endif  // DCA_LINALG_UTIL_CAST_CUDA_HPP
