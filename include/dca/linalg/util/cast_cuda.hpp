// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides the castCudaComplex utility and the CudaComplex typedef.

#ifndef DCA_LINALG_UTIL_CAST_CUDA_HPP
#define DCA_LINALG_UTIL_CAST_CUDA_HPP

#include <complex>
#include <cuComplex.h>

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

// Provides a templated typedef.
namespace details {
// dca::linalg::util::details::
template<typename Real>
struct ComplexContainer;
template<>
struct ComplexContainer<double>{
  using type = cuDoubleComplex;
};
template<>
struct ComplexContainer<float>{
  using type = cuFloatComplex;
};
}  // details
// dca::linalg::util::

template<typename Real>
using CudaComplex = typename details::ComplexContainer<Real>::type;

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_CAST_CUDA_HPP
