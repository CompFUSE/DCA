// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements laset_gpu.hpp.

#include "dca/linalg/lapack/multiply_diagonal_gpu.hpp"
#include <cassert>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace linalg {
namespace lapack {
namespace kernels {
// dca::linalg::lapack::kernels::

constexpr int multiply_diag_block_size_x = 32;
constexpr int multiply_diag_block_size_y = 32;

template <typename ScalarIn, typename ScalarOut>
__global__ void multiplyDiagonalLeft(int m, int n, const ScalarIn* d, int inc_d, const ScalarIn* a,
                                     int lda, ScalarOut* b, int ldb) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < m && j < n) {
    b[i + j * ldb] = d[i * inc_d] * a[i + j * lda];
  }
}

template <typename Type>
__global__ void multiplyDiagonalRight(int m, int n, const Type* a, int lda, const Type* d,
                                      int inc_d, Type* b, int ldb) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;

  if (i < m && j < n) {
    b[i + j * ldb] = d[j * inc_d] * a[i + j * lda];
  }
}

}  // namespace kernels
// dca::linalg::lapack::

template <typename ScalarIn, typename ScalarOut>
void multiplyDiagonalLeft_gpu(int m, int n, const ScalarIn* d, int inc_d, const ScalarIn* a,
                              int lda, ScalarOut* b, int ldb, int thread_id, int stream_id) {
  assert(lda >= m);
  assert(ldb >= m);

  if (m > 0 && n > 0) {
    checkErrorsCudaDebug();
    const int bl_x = dca::util::ceilDiv(m, kernels::multiply_diag_block_size_x);
    const int bl_y = dca::util::ceilDiv(n, kernels::multiply_diag_block_size_y);

    const dim3 threads(kernels::multiply_diag_block_size_x, kernels::multiply_diag_block_size_y);
    const dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::multiplyDiagonalLeft<ScalarIn, ScalarOut>
        <<<blocks, threads, 0, stream>>>(m, n, d, inc_d, a, lda, b, ldb);
    checkErrorsCudaDebug();
  }
}

template void multiplyDiagonalLeft_gpu<float, float>(int m, int n, const float* d, int inc_d,
                                                     const float* a, int lda, float* b, int ldb,
                                                     int thread_id, int stream_id);
template void multiplyDiagonalLeft_gpu<double, double>(int m, int n, const double* d, int inc_d,
                                                       const double* a, int lda, double* b, int ldb,
                                                       int thread_id, int stream_id);
template void multiplyDiagonalLeft_gpu<double, float>(int m, int n, const double* d, int inc_d,
                                                      const double* a, int lda, float* b, int ldb,
                                                      int thread_id, int stream_id);
template void multiplyDiagonalLeft_gpu<cuComplex, cuComplex>(int m, int n, const cuComplex* d,
                                                             int inc_d, const cuComplex* a, int lda,
                                                             cuComplex* b, int ldb, int thread_id,
                                                             int stream_id);
template void multiplyDiagonalLeft_gpu<cuDoubleComplex, cuDoubleComplex>(
    int m, int n, const cuDoubleComplex* d, int inc_d, const cuDoubleComplex* a, int lda,
    cuDoubleComplex* b, int ldb, int thread_id, int stream_id);

template <typename Type>
void multiplyDiagonalRight_gpu(int m, int n, const Type* a, int lda, const Type* d, int inc_d,
                               Type* b, int ldb, int thread_id, int stream_id) {
  assert(lda >= m);
  assert(ldb >= m);

  if (m > 0 && n > 0) {
    checkErrorsCudaDebug();
    const int bl_x = dca::util::ceilDiv(m, kernels::multiply_diag_block_size_x);
    const int bl_y = dca::util::ceilDiv(n, kernels::multiply_diag_block_size_y);

    const dim3 threads(kernels::multiply_diag_block_size_x, kernels::multiply_diag_block_size_y);
    const dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::multiplyDiagonalRight<<<blocks, threads, 0, stream>>>(m, n, a, lda, d, inc_d, b, ldb);
    checkErrorsCudaDebug();
  }
}
template void multiplyDiagonalRight_gpu(int m, int n, const float* a, int lda, const float* d,
                                        int inc_d, float* b, int ldb, int thread_id, int stream_id);
template void multiplyDiagonalRight_gpu(int m, int n, const double* a, int lda, const double* d,
                                        int inc_d, double* b, int ldb, int thread_id, int stream_id);
template void multiplyDiagonalRight_gpu(int m, int n, const cuComplex* a, int lda,
                                        const cuComplex* d, int inc_d, cuComplex* b, int ldb,
                                        int thread_id, int stream_id);
template void multiplyDiagonalRight_gpu(int m, int n, const cuDoubleComplex* a, int lda,
                                        const cuDoubleComplex* d, int inc_d, cuDoubleComplex* b,
                                        int ldb, int thread_id, int stream_id);

}  // namespace lapack
}  // namespace linalg
}  // namespace dca
