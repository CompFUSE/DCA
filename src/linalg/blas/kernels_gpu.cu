// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements kernels_gpu.hpp

#include "dca/linalg/blas/kernels_gpu.hpp"
#include <cassert>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace linalg {
namespace blas {
namespace kernels {
// dca::linalg::blas::kernels::

constexpr int copy_block_size_x = 32;
constexpr int copy_block_size_y = 8;
constexpr int scale_block_size_x = 32;
constexpr int scale_block_size_y = 8;

template <typename Type>
__global__ void copyRows(int row_size, int n_rows, const int* i_x, const Type* x, int ldx,
                         const int* i_y, Type* y, int ldy) {
  assert(blockDim.y == 1);
  assert(blockDim.z == 1);
  assert(blockIdx.z == 0);

  // Work on BlockDim.x rows and copyrows_block_size_y cols.
  int ind_i = threadIdx.x + blockIdx.x * blockDim.x;

  int js = blockIdx.y * copy_block_size_y;
  int je = min(row_size, (blockIdx.y + 1) * copy_block_size_y);

  if (ind_i < n_rows) {
    int iy = i_y[ind_i];
    int ix = i_x[ind_i];

    for (int j = js; j < je; ++j)
      y[iy + j * ldy] = x[ix + j * ldx];
  }
}

template <typename Type>
__global__ void copyCols(int col_size, int n_cols, const int* j_x, const Type* x, int ldx,
                         const int* j_y, Type* y, int ldy) {
  assert(blockDim.y == 1);
  assert(blockDim.z == 1);
  assert(blockIdx.z == 0);

  // Work on BlockDim.x rows and copyrows_block_size_y cols.
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  int ind_js = blockIdx.y * copy_block_size_y;
  int ind_je = min(n_cols, (blockIdx.y + 1) * copy_block_size_y);

  if (i < col_size) {
    for (int ind_j = ind_js; ind_j < ind_je; ++ind_j)
      y[i + j_y[ind_j] * ldy] = x[i + j_x[ind_j] * ldx];
  }
}

template <typename Type>
__global__ void scaleRows(int row_size, int n_rows, const int* i, const Type* alpha, Type* a,
                          int lda) {
  assert(blockDim.y == 1);
  assert(blockDim.z == 1);
  assert(blockIdx.z == 0);

  // Work on BlockDim.x rows and copyrows_block_size_y cols.
  int ind_i = threadIdx.x + blockIdx.x * blockDim.x;

  int js = blockIdx.y * scale_block_size_y;
  int je = min(row_size, (blockIdx.y + 1) * scale_block_size_y);

  if (ind_i < n_rows) {
    int ia = i[ind_i];

    for (int j = js; j < je; ++j)
      a[ia + j * lda] = a[ia + j * lda] * alpha[ind_i];
  }
}
}
// dca::linalg::blas::

template <typename Type>
void copyRows(int row_size, int n_rows, const int* i_x, const Type* x, int ldx, const int* i_y,
              Type* y, int ldy, int thread_id, int stream_id) {
  if (row_size > 0 && n_rows > 0) {
    checkErrorsCudaDebug();
    int bl_x = dca::util::ceilDiv(n_rows, kernels::copy_block_size_x);
    int bl_y = dca::util::ceilDiv(row_size, kernels::copy_block_size_y);

    dim3 threads(kernels::copy_block_size_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::copyRows<<<blocks, threads, 0, stream>>>(row_size, n_rows, i_x, x, ldx, i_y, y, ldy);
    checkErrorsCudaDebug();
  }
}
template void copyRows(int row_size, int n_rows, const int* i_x, const float* x, int ldx,
                       const int* i_y, float* y, int ldy, int thread_id, int stream_id);
template void copyRows(int row_size, int n_rows, const int* i_x, const double* x, int ldx,
                       const int* i_y, double* y, int ldy, int thread_id, int stream_id);
template void copyRows(int row_size, int n_rows, const int* i_x, const cuComplex* x, int ldx,
                       const int* i_y, cuComplex* y, int ldy, int thread_id, int stream_id);
template void copyRows(int row_size, int n_rows, const int* i_x, const cuDoubleComplex* x, int ldx,
                       const int* i_y, cuDoubleComplex* y, int ldy, int thread_id, int stream_id);

template <typename Type>
void copyCols(int col_size, int n_cols, const int* j_x, const Type* x, int ldx, const int* j_y,
              Type* y, int ldy, int thread_id, int stream_id) {
  if (col_size > 0 && n_cols > 0) {
    checkErrorsCudaDebug();
    int bl_x = dca::util::ceilDiv(col_size, kernels::copy_block_size_x);
    int bl_y = dca::util::ceilDiv(n_cols, kernels::copy_block_size_y);

    dim3 threads(kernels::copy_block_size_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::copyCols<<<blocks, threads, 0, stream>>>(col_size, n_cols, j_x, x, ldx, j_y, y, ldy);
    checkErrorsCudaDebug();
  }
}
template void copyCols(int col_size, int n_cols, const int* j_x, const float* x, int ldx,
                       const int* j_y, float* y, int ldy, int thread_id, int stream_id);
template void copyCols(int col_size, int n_cols, const int* j_x, const double* x, int ldx,
                       const int* j_y, double* y, int ldy, int thread_id, int stream_id);
template void copyCols(int col_size, int n_cols, const int* j_x, const cuComplex* x, int ldx,
                       const int* j_y, cuComplex* y, int ldy, int thread_id, int stream_id);
template void copyCols(int col_size, int n_cols, const int* j_x, const cuDoubleComplex* x, int ldx,
                       const int* j_y, cuDoubleComplex* y, int ldy, int thread_id, int stream_id);

template <typename Type>
void scaleRows(int row_size, int n_rows, const int* i, const Type* alpha, Type* a, int lda,
               int thread_id, int stream_id) {
  if (row_size > 0 && n_rows > 0) {
    checkErrorsCudaDebug();
    int bl_x = dca::util::ceilDiv(n_rows, kernels::scale_block_size_x);
    int bl_y = dca::util::ceilDiv(row_size, kernels::scale_block_size_y);

    dim3 threads(kernels::scale_block_size_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::scaleRows<<<blocks, threads, 0, stream>>>(row_size, n_rows, i, alpha, a, lda);
    checkErrorsCudaDebug();
  }
}
template void scaleRows(int row_size, int n_rows, const int* i, const float* alpha, float* a,
                        int lda, int thread_id, int stream_id);
template void scaleRows(int row_size, int n_rows, const int* i, const double* alpha, double* a,
                        int lda, int thread_id, int stream_id);
template void scaleRows(int row_size, int n_rows, const int* i, const cuComplex* alpha,
                        cuComplex* a, int lda, int thread_id, int stream_id);
template void scaleRows(int row_size, int n_rows, const int* i, const cuDoubleComplex* alpha,
                        cuDoubleComplex* a, int lda, int thread_id, int stream_id);
}  // blas
}  // linalg
}  // dca
