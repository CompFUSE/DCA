// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides device-generic wrappers for some of the LAPACK routines.

#ifndef DCA_LINALG_LAPACK_USE_DEVICE_HPP
#define DCA_LINALG_LAPACK_USE_DEVICE_HPP

#include "dca/linalg/device_type.hpp"

#include "dca/linalg/lapack/lapack.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/lapack/laset_gpu.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/lapack/multiply_diagonal_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

template <DeviceType DeviceName>
struct UseDevice;

template <>
struct UseDevice<CPU> {
  // Auxiliary routines
  template <typename ScalarType>
  inline static void laset(int m, int n, ScalarType offdiag, ScalarType diag, ScalarType* a,
                           int lda, int /*thread_id*/, int /*stream_id*/) {
    lapack::laset("A", m, n, offdiag, diag, a, lda);
  }

  // Computational routines
  template <typename ScalarType>
  inline static void getrf(int m, int n, ScalarType* a, int lda, int* ipiv) {
    lapack::getrf(m, n, a, lda, ipiv);
  }

  template <typename ScalarType>
  inline static void getri(int n, ScalarType* a, int lda, int* ipiv, ScalarType* work, int lwork) {
    lapack::getri(n, a, lda, ipiv, work, lwork);
  }

  // Custom routines

  // Performs the matrix-matrix multiplication b <- d * a,
  // where d is a diagonal matrix, which diagonal elements are given by d[0], d[inc_d],
  // d[2*inc_d]...
  // Out: b
  // Preconditions: lda >= m, ldb >= m.
  template <typename ScalarIn, typename ScalarOut>
  static void multiplyDiagonalLeft(int m, int n, const ScalarIn* d, int inc_d,
                                   const ScalarIn* a, int lda, ScalarOut* b, int ldb,
                                   int /*thread_id*/, int /*stream_id*/) {
    assert(lda >= m);
    assert(ldb >= m);

    for (int j = 0; j < n; ++j) {
      const ScalarIn* aj = a + j * lda;
      ScalarOut* bj = b + j * ldb;

      for (int i = 0; i < m; ++i)
        bj[i] = aj[i] * d[i * inc_d];
    }
  }

  // Performs the matrix-matrix multiplication b <- a * d,
  // where d is a diagonal matrix, which diagonal elements are given by d[0], d[inc_d],
  // d[2*inc_d]...
  // Out: b
  // Preconditions: lda >= m, ldb >= m.
  template <typename ScalarType>
  static void multiplyDiagonalRight(int m, int n, const ScalarType* a, int lda, const ScalarType* d,
                                    int inc_d, ScalarType* b, int ldb, int /*thread_id*/,
                                    int /*stream_id*/) {
    assert(lda >= m);
    assert(ldb >= m);

    for (int j = 0; j < n; ++j) {
      ScalarType dj = d[j * inc_d];
      const ScalarType* aj = a + j * lda;
      ScalarType* bj = b + j * ldb;

      for (int i = 0; i < m; ++i)
        bj[i] = aj[i] * dj;
    }
  }
};

#ifdef DCA_HAVE_CUDA
template <>
struct UseDevice<GPU> {
  // Auxiliary routines
  template <typename ScalarType>
  inline static void laset(int m, int n, ScalarType offdiag, ScalarType diag, ScalarType* a,
                           int lda, int thread_id, int stream_id) {
    lapack::laset_gpu(m, n, offdiag, diag, a, lda, thread_id, stream_id);
  }

  // Computational routines
  template <typename ScalarType>
  inline static void getrf(int m, int n, ScalarType* a, int lda, int* ipiv) {
    magma::getrf_gpu(m, n, a, lda, ipiv);
  }

  template <typename ScalarType>
  inline static void getri(int n, ScalarType* a, int lda, int* ipiv, ScalarType* work, int lwork) {
    magma::getri_gpu(n, a, lda, ipiv, work, lwork);
  }

  template <typename ScalarIn, typename ScalarOut>
  inline static void multiplyDiagonalLeft(int m, int n, const ScalarIn* d, int inc_d,
                                          const ScalarIn* a, int lda, ScalarOut* b, int ldb,
                                          int thread_id, int stream_id) {
    lapack::multiplyDiagonalLeft_gpu(m, n, d, inc_d, a, lda, b, ldb, thread_id, stream_id);
  }

  template <typename ScalarType>
  inline static void multiplyDiagonalRight(int m, int n, const ScalarType* a, int lda,
                                           const ScalarType* d, int inc_d, ScalarType* b, int ldb,
                                           int thread_id, int stream_id) {
    lapack::multiplyDiagonalRight_gpu(m, n, a, lda, d, inc_d, b, ldb, thread_id, stream_id);
  }
};
#endif  // DCA_HAVE_CUDA

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_USE_DEVICE_HPP
