// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides device-generic wrappers for some of the BLAS routines.

#ifndef DCA_LINALG_BLAS_USE_DEVICE_HPP
#define DCA_LINALG_BLAS_USE_DEVICE_HPP

#include "dca/linalg/device_type.hpp"

#include "dca/linalg/blas/blas1.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/blas/blas3.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/blas/cublas1.hpp"
#include "dca/linalg/blas/cublas3.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace blas {
// dca::linalg::blas::

template <DeviceType DeviceName>
struct UseDevice;

template <>
struct UseDevice<CPU> {
  // Level 1
  template <typename ScalarType>
  inline static void axpy(int n, ScalarType alpha, ScalarType* x, int incx, ScalarType* y, int incy,
                          int /*thread_id*/, int /*stream_id*/) {
    blas::axpy(n, alpha, x, incx, y, incy);
  }

  // Level 3
  template <typename ScalarType>
  inline static void gemm(const char* transa, const char* transb, int m, int n, int k,
                          ScalarType alpha, const ScalarType* a, int lda, const ScalarType* b,
                          int ldb, ScalarType beta, ScalarType* c, int ldc, int /*thread_id*/,
                          int /*stream_id*/) {
    blas::gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  }

  template <typename ScalarType>
  inline static void trsm(const char* side, const char* uplo, const char* transa, const char* diag,
                          int m, int n, ScalarType alpha, const ScalarType* a, int lda,
                          ScalarType* b, int ldb, int /*thread_id*/, int /*stream_id*/) {
    blas::trsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  }
};

#ifdef DCA_HAVE_CUDA
using LIN_ALG::get_thread_handle;

template <>
struct UseDevice<GPU> {
  // Level 1
  template <typename ScalarType>
  inline static void axpy(int n, ScalarType alpha, const ScalarType* x, int incx, ScalarType* y,
                          int incy, int thread_id, int /*stream_id*/) {
    cublas::axpy(get_thread_handle(thread_id), n, alpha, x, incx, y, incy);
  }

  // Level 3
  template <typename ScalarType>
  inline static void gemm(const char* transa, const char* transb, int m, int n, int k,
                          ScalarType alpha, const ScalarType* a, int lda, const ScalarType* b,
                          int ldb, ScalarType beta, ScalarType* c, int ldc, int thread_id,
                          int /*stream_id*/) {
    cublas::gemm(get_thread_handle(thread_id), transa, transb, m, n, k, alpha, a, lda, b, ldb, beta,
                 c, ldc);
  }

  template <typename ScalarType>
  inline static void trsm(const char* side, const char* uplo, const char* transa, const char* diag,
                          int m, int n, ScalarType alpha, const ScalarType* a, int lda,
                          ScalarType* b, int ldb, int thread_id, int /*stream_id*/) {
    cublas::trsm(get_thread_handle(thread_id), side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  }
};
#endif  // DCA_HAVE_CUDA

}  // blas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_USE_DEVICE_HPP
