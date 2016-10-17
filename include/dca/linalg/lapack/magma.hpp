// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the C++ wrappers for some of the MAGMA functions.

#ifndef DCA_LINALG_LAPACK_MAGMA_HPP
#define DCA_LINALG_LAPACK_MAGMA_HPP

#include <complex>
#include <magma.h>

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/error_cuda.hpp"

// C++ wrappers
namespace dca {
namespace linalg {
namespace magma {
// dca::linalg::magma::

inline void getrf_gpu(int m, int n, float* a, int lda, int* ipiv, int* info) {
  checkErrorsCudaDebug();
  magma_sgetrf_gpu(m, n, a, lda, ipiv, info);
  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, double* a, int lda, int* ipiv, int* info) {
  checkErrorsCudaDebug();
  magma_dgetrf_gpu(m, n, a, lda, ipiv, info);
  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, std::complex<float>* a, int lda, int* ipiv, int* info) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  magma_cgetrf_gpu(m, n, cu_a, lda, ipiv, info);
  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, std::complex<double>* a, int lda, int* ipiv, int* info) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  magma_zgetrf_gpu(m, n, cu_a, lda, ipiv, info);
  checkErrorsCudaDebug();
}

inline void getri_gpu(int n, float* a, int lda, int* ipiv, float* work, int lwork, int* info) {
  checkErrorsCudaDebug();
  magma_sgetri_gpu(n, a, lda, ipiv, work, lwork, info);
  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, double* a, int lda, int* ipiv, double* work, int lwork, int* info) {
  checkErrorsCudaDebug();
  magma_dgetri_gpu(n, a, lda, ipiv, work, lwork, info);
  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, std::complex<float>* a, int lda, int* ipiv, std::complex<float>* work,
                      int lwork, int* info) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  auto cu_work = util::castCudaComplex(work);
  magma_cgetri_gpu(n, cu_a, lda, ipiv, cu_work, lwork, info);
  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, std::complex<double>* a, int lda, int* ipiv,
                      std::complex<double>* work, int lwork, int* info) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  auto cu_work = util::castCudaComplex(work);
  magma_zgetri_gpu(n, cu_a, lda, ipiv, cu_work, lwork, info);
  checkErrorsCudaDebug();
}

template <typename Type>
int get_getri_nb(int n);
template <>
inline int get_getri_nb<float>(int n) {
  return magma_get_sgetri_nb(n);
}
template <>
inline int get_getri_nb<double>(int n) {
  return magma_get_dgetri_nb(n);
}
template <>
inline int get_getri_nb<std::complex<float>>(int n) {
  return magma_get_cgetri_nb(n);
}
template <>
inline int get_getri_nb<std::complex<double>>(int n) {
  return magma_get_zgetri_nb(n);
}

}  // magma
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_MAGMA_HPP
