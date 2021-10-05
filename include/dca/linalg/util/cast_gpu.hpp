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

#if defined(DCA_HAVE_CUDA)
#include <cuComplex.h>
#elif defined(DCA_HAVE_HIP)
#include <hip/hip_complex.h>
#include "dca/util/cuda2hip.h"
#endif

#include <magma_v2.h>
namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

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
  using type = cuFloatComplex;
};
}  // namespace details
// dca::linalg::util::

template <typename Real>
using CudaComplex = typename details::ComplexContainer<Real>::type;

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_CAST_CUDA_HPP
