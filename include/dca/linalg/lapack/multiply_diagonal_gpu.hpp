// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the GPU implementation of the laset function.

#ifndef DCA_LINALG_LAPACK_MULTIPLY_DIAGONAL_GPU_HPP
#define DCA_LINALG_LAPACK_MULTIPLY_DIAGONAL_GPU_HPP

#include <complex>
#include <cuComplex.h>
#include "dca/linalg/util/cast_cuda.hpp"

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

// Performs the matrix-matrix multiplication b <- d * a,
// where d is a diagonal matrix, which diagonal elements are given by d[0], d[inc_d], d[2*inc_d]...
// Out: b
// Preconditions: lda >= m, ldb >= m.
// Type can be float, double, cuComplex, cuDoubleComplex, std::complex<float>, std::complex<double>.
template <typename ScalarIn, typename ScalarOut>
void multiplyDiagonalLeft_gpu(int m, int n, const ScalarIn* d, int inc_d, const ScalarIn* a,
                              int lda, ScalarOut* b, int ldb, int thread_id, int stream_id);
template <typename Type>
inline void multiplyDiagonalLeft_gpu(int m, int n, const std::complex<Type>* d, int inc_d,
                                     const std::complex<Type>* a, int lda, std::complex<Type>* b,
                                     int ldb, int thread_id, int stream_id) {
  auto cu_a = linalg::util::castCudaComplex(a);
  auto cu_d = linalg::util::castCudaComplex(d);
  auto cu_b = linalg::util::castCudaComplex(b);
  multiplyDiagonalLeft_gpu(m, n, cu_d, inc_d, cu_a, lda, cu_b, ldb, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication b <- d * a,
// where d is a diagonal matrix, which diagonal elements are given by d[0], d[inc_d], d[2*inc_d]...
// Out: b
// Preconditions: lda >= m, ldb >= m.
// Type can be float, double, cuComplex, cuDoubleComplex, std::complex<float>, std::complex<double>.
template <typename Type>
void multiplyDiagonalRight_gpu(int m, int n, const Type* a, int lda, const Type* d, int inc_d,
                               Type* b, int ldb, int thread_id, int stream_id);
template <typename Type>
inline void multiplyDiagonalRight_gpu(int m, int n, const std::complex<Type>* a, int lda,
                                      const std::complex<Type>* d, int inc_d, std::complex<Type>* b,
                                      int ldb, int thread_id, int stream_id) {
  auto cu_a = linalg::util::castCudaComplex(a);
  auto cu_d = linalg::util::castCudaComplex(d);
  auto cu_b = linalg::util::castCudaComplex(b);
  multiplyDiagonalRight_gpu(m, n, cu_a, lda, cu_d, inc_d, cu_b, ldb, thread_id, stream_id);
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_MULTIPLY_DIAGONAL_GPU_HPP
