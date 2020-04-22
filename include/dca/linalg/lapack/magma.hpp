// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the C++ wrappers for some of the MAGMA functions.

#ifndef DCA_LINALG_LAPACK_MAGMA_HPP
#define DCA_LINALG_LAPACK_MAGMA_HPP

#include <cassert>
#include <complex>
#include <magma.h>

#include "dca/linalg/lapack/lapack.hpp"
#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/error_cuda.hpp"

// C++ wrappers
namespace dca {
namespace linalg {
namespace magma {
// dca::linalg::magma::

inline void getrf_gpu(int m, int n, float* a, int lda, int* ipiv) {
  checkErrorsCudaDebug();

  int info = 0;
  magma_sgetrf_gpu(m, n, a, lda, ipiv, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, double* a, int lda, int* ipiv) {
  checkErrorsCudaDebug();

  int info = 0;
  magma_dgetrf_gpu(m, n, a, lda, ipiv, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, std::complex<float>* a, int lda, int* ipiv) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);

  int info = 0;
  magma_cgetrf_gpu(m, n, cu_a, lda, ipiv, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getrf_gpu(int m, int n, std::complex<double>* a, int lda, int* ipiv) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);

  int info = 0;
  magma_zgetrf_gpu(m, n, cu_a, lda, ipiv, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}

inline void getri_gpu(int n, float* a, int lda, int* ipiv, float* work, int lwork) {
  checkErrorsCudaDebug();

  int info = 0;
  magma_sgetri_gpu(n, a, lda, ipiv, work, lwork, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, double* a, int lda, int* ipiv, double* work, int lwork) {
  checkErrorsCudaDebug();

  int info = 0;
  magma_dgetri_gpu(n, a, lda, ipiv, work, lwork, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, std::complex<float>* a, int lda, int* ipiv, std::complex<float>* work,
                      int lwork) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  auto cu_work = util::castCudaComplex(work);

  int info = 0;
  magma_cgetri_gpu(n, cu_a, lda, ipiv, cu_work, lwork, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}
inline void getri_gpu(int n, std::complex<double>* a, int lda, int* ipiv,
                      std::complex<double>* work, int lwork) {
  checkErrorsCudaDebug();
  auto cu_a = util::castCudaComplex(a);
  auto cu_work = util::castCudaComplex(work);

  int info = 0;
  magma_zgetri_gpu(n, cu_a, lda, ipiv, cu_work, lwork, &info);
  checkLapackInfo(info);

  checkErrorsCudaDebug();
}

inline magma_trans_t toMagmaTrans(const char x) {
  switch (x) {
    case 'N':
      return MagmaNoTrans;
    case 'T':
      return MagmaTrans;
    case 'C':
      return MagmaConjTrans;
    default:
      throw(std::logic_error(__FUNCTION__));
  }
}

inline void magmablas_gemm_vbatched(const char transa, const char transb, int* m, int* n, int* k,
                                    const float alpha, const float* const* a, int* lda,
                                    const float* const* b, int* ldb, const float beta, float** c,
                                    int* ldc, const int batch_count, magma_queue_t queue) {
  magmablas_sgemm_vbatched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha, a, lda, b,
                           ldb, beta, c, ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_vbatched(const char transa, const char transb, int* m, int* n, int* k,
                                    const double alpha, const double* const* a, int* lda,
                                    const double* const* b, int* ldb, const double beta, double** c,
                                    int* ldc, const int batch_count, const magma_queue_t queue) {
  magmablas_dgemm_vbatched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha, a, lda, b,
                           ldb, beta, c, ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_vbatched(const char transa, const char transb, int* m, int* n, int* k,
                                    const std::complex<float> alpha,
                                    const std::complex<float>* const* a, int* lda,
                                    const std::complex<float>* const* b, int* ldb,
                                    const std::complex<float> beta, std::complex<float>** c,
                                    int* ldc, const int batch_count, const magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_cgemm_vbatched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k,
                           *castCudaComplex(alpha), castCudaComplex(a), lda, castCudaComplex(b),
                           ldb, *castCudaComplex(beta), castCudaComplex(c), ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_vbatched(const char transa, const char transb, int* m, int* n, int* k,
                                    const std::complex<double> alpha,
                                    const std::complex<double>* const* a, int* lda,
                                    const std::complex<double>* const* b, int* ldb,
                                    const std::complex<double> beta, std::complex<double>** c,
                                    int* ldc, const int batch_count, const magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_zgemm_vbatched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k,
                           *castCudaComplex(alpha), castCudaComplex(a), lda, castCudaComplex(b),
                           ldb, *castCudaComplex(beta), castCudaComplex(c), ldc, batch_count, queue);
  checkErrorsCudaDebug();
}

inline void magmablas_gemm_vbatched_max_nocheck(const char transa, const char transb, int* m,
                                                int* n, int* k, const float alpha,
                                                const float* const* a, int* lda,
                                                const float* const* b, int* ldb, const float beta,
                                                float** c, int* ldc, const int batch_count,
                                                const int m_max, const int n_max, const int k_max,
                                                magma_queue_t queue) {
  magmablas_sgemm_vbatched_max_nocheck(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha,
                                       a, lda, b, ldb, beta, c, ldc, batch_count, m_max, n_max,
                                       k_max, queue);
  checkErrorsCudaDebug();
}

inline void magmablas_gemm_vbatched_max_nocheck(const char transa, const char transb, int* m,
                                                int* n, int* k, const double alpha,
                                                const double* const* a, int* lda,
                                                const double* const* b, int* ldb, const double beta,
                                                double** c, int* ldc, const int batch_count,
                                                const int m_max, const int n_max, const int k_max,
                                                magma_queue_t queue) {
  magmablas_dgemm_vbatched_max_nocheck(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha,
                                       a, lda, b, ldb, beta, c, ldc, batch_count, m_max, n_max,
                                       k_max, queue);
  checkErrorsCudaDebug();
}

inline void magmablas_gemm_vbatched_max_nocheck(
    const char transa, const char transb, int* m, int* n, int* k, const std::complex<float> alpha,
    const std::complex<float>* const* a, int* lda, const std::complex<float>* const* b, int* ldb,
    const std::complex<float> beta, std::complex<float>** c, int* ldc, const int batch_count,
    const int m_max, const int n_max, const int k_max, magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_cgemm_vbatched_max_nocheck(
      toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, *castCudaComplex(alpha),
      castCudaComplex(a), lda, castCudaComplex(b), ldb, *castCudaComplex(beta), castCudaComplex(c),
      ldc, batch_count, m_max, n_max, k_max, queue);
  checkErrorsCudaDebug();
}

inline void magmablas_gemm_vbatched_max_nocheck(
    const char transa, const char transb, int* m, int* n, int* k, const std::complex<double> alpha,
    const std::complex<double>* const* a, int* lda, const std::complex<double>* const* b, int* ldb,
    const std::complex<double> beta, std::complex<double>** c, int* ldc, const int batch_count,
    const int m_max, const int n_max, const int k_max, magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_zgemm_vbatched_max_nocheck(
      toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, *castCudaComplex(alpha),
      castCudaComplex(a), lda, castCudaComplex(b), ldb, *castCudaComplex(beta), castCudaComplex(c),
      ldc, batch_count, m_max, n_max, k_max, queue);
  checkErrorsCudaDebug();
}

inline void magmablas_gemm_batched(const char transa, const char transb, const int m, const int n,
                                   const int k, const float alpha, const float* const* a,
                                   const int lda, const float* const* b, const int ldb,
                                   const float beta, float** c, const int ldc,
                                   const int batch_count, magma_queue_t queue) {
  magmablas_sgemm_batched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha, a, lda, b,
                          ldb, beta, c, ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_batched(const char transa, const char transb, const int m, const int n,
                                   const int k, const double alpha, const double* const* a,
                                   const int lda, const double* const* b, const int ldb,
                                   const double beta, double** c, const int ldc,
                                   const int batch_count, const magma_queue_t queue) {
  magmablas_dgemm_batched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k, alpha, a, lda, b,
                          ldb, beta, c, ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_batched(const char transa, const char transb, const int m, const int n,
                                   const int k, const std::complex<float> alpha,
                                   const std::complex<float>* const* a, const int lda,
                                   const std::complex<float>* const* b, const int ldb,
                                   const std::complex<float> beta, std::complex<float>** c,
                                   const int ldc, const int batch_count, const magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_cgemm_batched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k,
                          *castCudaComplex(alpha), castCudaComplex(a), lda, castCudaComplex(b), ldb,
                          *castCudaComplex(beta), castCudaComplex(c), ldc, batch_count, queue);
  checkErrorsCudaDebug();
}
inline void magmablas_gemm_batched(const char transa, const char transb, const int m, const int n,
                                   const int k, const std::complex<double> alpha,
                                   const std::complex<double>* const* a, const int lda,
                                   const std::complex<double>* const* b, const int ldb,
                                   const std::complex<double> beta, std::complex<double>** c,
                                   const int ldc, const int batch_count, const magma_queue_t queue) {
  using util::castCudaComplex;
  magmablas_zgemm_batched(toMagmaTrans(transa), toMagmaTrans(transb), m, n, k,
                          *castCudaComplex(alpha), castCudaComplex(a), lda, castCudaComplex(b), ldb,
                          *castCudaComplex(beta), castCudaComplex(c), ldc, batch_count, queue);
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

}  // namespace magma
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_LAPACK_MAGMA_HPP
