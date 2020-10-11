// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides the castCudaComplex utility and the CudaComplex typedef.

#ifndef DCA_LINALG_UTIL_CAST_CUDA_HPP
#define DCA_LINALG_UTIL_CAST_CUDA_HPP

#include <complex>
#include <cuComplex.h>

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

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

template <class T>
struct CudaConvert {
  using type = T;
};

template <class T>
struct CudaConvert<std::complex<T>> {
  using type = typename ComplexContainer<T>::type;
};
}  // namespace details
// dca::linalg::util::

template <typename Real>
using CudaComplex = typename details::ComplexContainer<Real>::type;

template <typename T>
using CudaConvert = typename details::CudaConvert<T>::type;

template <class T>
__device__ __host__ static auto& castCuda(T& x) {
  return reinterpret_cast<CudaConvert<T>&>(x);
}

template <class T>
__device__ __host__ static const auto& castCuda(const T& x) {
  return reinterpret_cast<const CudaConvert<T>&>(x);
}

template <class T>
__device__ __host__ static auto* castCuda(T* x) {
  return reinterpret_cast<CudaConvert<T>*>(x);
}

template <class T>
__device__ __host__ static const auto* castCuda(const T* x) {
  return reinterpret_cast<const CudaConvert<T>*>(x);
}

template <class T>
__device__ __host__ static auto* castCuda(T** x) {
  return reinterpret_cast<CudaConvert<T>**>(x);
}

template <class T>
__device__ __host__ static auto* castCuda(T const* const* x) {
  return reinterpret_cast<CudaConvert<T> const* const*>(x);
}

template <class T>
__device__ __host__ static auto& castCuda(MatrixView<T, GPU>& x) {
  return reinterpret_cast<MatrixView<CudaConvert<T>, GPU>&>(x);
}

template <class T>
__device__ __host__ static const auto& castCuda(const MatrixView<T, GPU>& x) {
  return reinterpret_cast<const MatrixView<CudaConvert<T>, GPU>&>(x);
}

}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_CAST_CUDA_HPP
