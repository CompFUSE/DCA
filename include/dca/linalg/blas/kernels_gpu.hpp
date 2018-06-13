// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides some BLAS-like functions.

#ifndef DCA_LINALG_BLAS_KERNELS_GPU_HPP
#define DCA_LINALG_BLAS_KERNELS_GPU_HPP

#include <complex>
#include <cuComplex.h>
#include "dca/linalg/util/cast_cuda.hpp"

namespace dca {
namespace linalg {
namespace blas {
// dca::linalg::blas::

// Type can be float, double, cuComplex, cuDoubleComplex, std::complex<float>, std::complex<double>.
template <typename Type>
void copyRows(int row_size, int n_rows, const int* i_x, const Type* x, int ldx, const int* i_y,
              Type* y, int ldy, int thread_id, int stream_id);
template <typename Type>
inline void copyRows(int row_size, int n_rows, const int* i_x, const std::complex<Type>* x, int ldx,
                     const int* i_y, std::complex<Type>* y, int ldy, int thread_id, int stream_id) {
  auto cu_x = util::castCudaComplex(x);
  auto cu_y = util::castCudaComplex(y);
  copyRows(row_size, n_rows, i_x, cu_x, ldx, i_y, cu_y, ldy, thread_id, stream_id);
}

template <typename Type>
void copyCols(int col_size, int n_cols, const int* j_x, const Type* x, int ldx, const int* j_y,
              Type* y, int ldy, int thread_id, int stream_id);
template <typename Type>
inline void copyCols(int col_size, int n_cols, const int* j_x, const std::complex<Type>* x, int ldx,
                     const int* j_y, std::complex<Type>* y, int ldy, int thread_id, int stream_id) {
  auto cu_x = util::castCudaComplex(x);
  auto cu_y = util::castCudaComplex(y);
  copyCols(col_size, n_cols, j_x, cu_x, ldx, j_y, cu_y, ldy, thread_id, stream_id);
}

template <typename Type>
void moveLeft(int m, int n, Type* a, int lda);
template <typename Type>
inline void moveLeft(int m, int n, std::complex<Type>* a, int lda) {
  auto cu_a = util::castCudaComplex(a);
  moveLeft(m, n, cu_a, lda);
}

template <typename Type>
void moveUp(int m, int n, Type* a, int lda);
template <typename Type>
inline void moveUp(int m, int n, std::complex<Type>* a, int lda) {
  auto cu_a = util::castCudaComplex(a);
  moveUp(m, n, cu_a, lda);
}

template <typename Type>
void scaleRows(int row_size, int n_rows, const int* i, const Type* alpha, Type* a, int lda,
               int thread_id, int stream_id);
template <typename Type>
inline void scaleRows(int row_size, int n_rows, const int* i, const std::complex<Type>* alpha,
                      std::complex<Type>* a, int lda, int thread_id, int stream_id) {
  auto cu_alpha = util::castCudaComplex(alpha);
  auto cu_a = util::castCudaComplex(a);
  scaleRows(row_size, n_rows, i, cu_alpha, cu_a, lda, thread_id, stream_id);
}

template <typename Type>
void swapRows(int row_size, int n_rows, const int* i_1, const int* i_2, Type* a, int lda,
              int thread_id, int stream_id);
template <typename Type>
inline void swapRows(int row_size, int n_rows, const int* i_1, const int* i_2,
                     std::complex<Type>* a, int lda, int thread_id, int stream_id) {
  auto cu_a = util::castCudaComplex(a);
  swapRows(row_size, n_rows, i_1, i_2, cu_a, lda, thread_id, stream_id);
}

template <typename Type>
void swapCols(int col_size, int n_cols, const int* j_1, const int* j_2, Type* a, int lda,
              int thread_id, int stream_id);
template <typename Type>
inline void swapCols(int col_size, int n_cols, const int* j_1, const int* j_2,
                     std::complex<Type>* a, int lda, int thread_id, int stream_id) {
  auto cu_a = util::castCudaComplex(a);
  swapCols(col_size, n_cols, j_1, j_2, cu_a, lda, thread_id, stream_id);
}
}  // blas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_KERNELS_GPU_HPP
