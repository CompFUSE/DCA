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
#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/util/type_help.hpp"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/cast_gpu.hpp"

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
void copyRows(int row_size, int n_rows, const int* i_x, const Type* x, int ldx, Type* y, int ldy,
              int thread_id, int stream_id);

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
void copyCols(int col_size, int n_cols, const int* j_x, const Type* x, int ldx, Type* y, int ldy,
              int thread_id, int stream_id);

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

extern template void copyCols(int col_size, int n_cols, const int* j_x, const float* x, int ldx,
                       const int* j_y, float* y, int ldy, int thread_id, int stream_id);
extern template void copyCols(int col_size, int n_cols, const int* j_x, const double* x, int ldx,
                       const int* j_y, double* y, int ldy, int thread_id, int stream_id);
extern template void copyCols(int col_size, int n_cols, const int* j_x, const cuComplex* x, int ldx,
                       const int* j_y, cuComplex* y, int ldy, int thread_id, int stream_id);
extern template void copyCols(int col_size, int n_cols, const int* j_x, const cuDoubleComplex* x, int ldx,
                       const int* j_y, cuDoubleComplex* y, int ldy, int thread_id, int stream_id);
extern template void copyCols(int col_size, int n_cols, const int* j_x, const double* x, int ldx,
                       const int* j_y, double* y, int ldy, int thread_id, int stream_id);
extern template void copyCols(int col_size, int n_cols, const int* j_x, const float* x, int ldx,
                       const int* j_y, float* y, int ldy, int thread_id, int stream_id);

extern template void copyCols(int col_size, int n_cols, const int* j_x, const cuComplex* x, int ldx,
                       const int* j_y, cuComplex* y, int ldy, int thread_id, int stream_id);

extern template void copyCols(int, int, const int*, const float*, int, float*, int, int, int);
extern template void copyCols(int, int, const int*, const double*, int, double*, int, int, int);

extern template void copyRows(int row_size, int n_rows, const int* i_x, const float* x, int ldx,
                       const int* i_y, float* y, int ldy, int thread_id, int stream_id);
extern template void copyRows(int row_size, int n_rows, const int* i_x, const double* x, int ldx,
                       const int* i_y, double* y, int ldy, int thread_id, int stream_id);
extern template void copyRows(int row_size, int n_rows, const int* i_x, const cuComplex* x, int ldx,
                       const int* i_y, cuComplex* y, int ldy, int thread_id, int stream_id);
extern template void copyRows(int row_size, int n_rows, const int* i_x, const cuDoubleComplex* x, int ldx,
                       const int* i_y, cuDoubleComplex* y, int ldy, int thread_id, int stream_id);
extern template void moveLeft(int m, int n, float* a, int lda);
extern template void moveLeft(int m, int n, double* a, int lda);
extern template void moveLeft(int m, int n, cuComplex* a, int lda);
extern template void moveLeft(int m, int n, cuDoubleComplex* a, int lda);
extern template void moveUp(int m, int n, float* a, int lda);
extern template void moveUp(int m, int n, double* a, int lda);
extern template void moveUp(int m, int n, cuComplex* a, int lda);
extern template void moveUp(int m, int n, cuDoubleComplex* a, int lda);
extern template void swapCols(int col_size, int n_cols, const int* j1, const int* j2, float* a, int lda,
                       int thread_id, int stream_id);
extern template void swapCols(int col_size, int n_cols, const int* j1, const int* j2, double* a, int lda,
                       int thread_id, int stream_id);
extern template void swapCols(int col_size, int n_cols, const int* j1, const int* j2, cuComplex* a,
                       int lda, int thread_id, int stream_id);
extern template void swapCols(int col_size, int n_cols, const int* j1, const int* j2, cuDoubleComplex* a,
                       int lda, int thread_id, int stream_id);
extern template void scaleRows(int row_size, int n_rows, const int* i, const float* alpha, float* a,
                        int lda, int thread_id, int stream_id);
extern template void scaleRows(int row_size, int n_rows, const int* i, const double* alpha, double* a,
                        int lda, int thread_id, int stream_id);
extern template void scaleRows(int row_size, int n_rows, const int* i, const cuComplex* alpha,
                        cuComplex* a, int lda, int thread_id, int stream_id);
extern template void scaleRows(int row_size, int n_rows, const int* i, const cuDoubleComplex* alpha,
                        cuDoubleComplex* a, int lda, int thread_id, int stream_id);
extern template void swapRows(int row_size, int n_rows, const int* i1, const int* i2, float* a, int lda,
                       int thread_id, int stream_id);
extern template void swapRows(int row_size, int n_rows, const int* i1, const int* i2, double* a, int lda,
                       int thread_id, int stream_id);
extern template void swapRows(int row_size, int n_rows, const int* i1, const int* i2, cuComplex* a,
                       int lda, int thread_id, int stream_id);
extern template void swapRows(int row_size, int n_rows, const int* i1, const int* i2, cuDoubleComplex* a,
                       int lda, int thread_id, int stream_id);



  
}  // namespace blas
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_BLAS_KERNELS_GPU_HPP
