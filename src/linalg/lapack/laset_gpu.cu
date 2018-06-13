// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements laset_gpu.hpp.

#include "dca/linalg/lapack/laset_gpu.hpp"
#include <cassert>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace linalg {
namespace lapack {
namespace kernels {
// dca::linalg::lapack::kernels::

constexpr int laset_block_size = 32;

template <typename Type>
__device__ void laset_set_diag(int m, int n, Type offdiag, Type diag, Type* a, int lda) {
  // Work on a tile of size (blockDim.x x blockDim.x).
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i < m) {
    int js = blockIdx.y * blockDim.x;
    int je = min(n, (blockIdx.y + 1) * blockDim.x);

    for (int j = js; j < je; ++j)
      if (i == j)
        a[i + j * lda] = diag;
      else
        a[i + j * lda] = offdiag;
  }
}

template <typename Type>
__device__ void laset_set_offdiag(int m, int n, Type offdiag, Type* a, int lda) {
  // Work on a tile of size (blockDim.x x blockDim.x).
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i < m) {
    int js = blockIdx.y * blockDim.x;
    int je = min(n, (blockIdx.y + 1) * blockDim.x);

    for (int j = js; j < je; ++j)
      a[i + j * lda] = offdiag;
  }
}

template <typename Type>
__global__ void laset(int m, int n, Type offdiag, Type diag, Type* a, int lda) {
  assert(blockDim.y == 1);
  assert(blockDim.z == 1);
  assert(blockIdx.z == 0);

  if (blockIdx.x == blockIdx.y) {
    laset_set_diag(m, n, offdiag, diag, a, lda);
  }
  else {
    laset_set_offdiag(m, n, offdiag, a, lda);
  }
}
}  // kernels
// dca::linalg::lapack::

// Sets the diagonal elements of the matrix to diag, and the off diagonal elements to offdiag.
// Preconditions: lda >= m.
// Type can be float, double, cuComplex, cuDoubleComplex.
template <typename Type>
void laset_gpu(int m, int n, Type offdiag, Type diag, Type* a, int lda, int thread_id, int stream_id) {
  assert(lda >= m);

  if (m > 0 && n > 0) {
    checkErrorsCudaDebug();
    int bl_x = dca::util::ceilDiv(m, kernels::laset_block_size);
    int bl_y = dca::util::ceilDiv(n, kernels::laset_block_size);

    dim3 threads(kernels::laset_block_size);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream = dca::linalg::util::getStream(thread_id, stream_id);

    kernels::laset<<<blocks, threads, 0, stream>>>(m, n, offdiag, diag, a, lda);
    checkErrorsCudaDebug();
  }
}

template void laset_gpu(int m, int n, float offdiag, float diag, float* a, int lda, int thread_id,
                        int stream_id);
template void laset_gpu(int m, int n, double offdiag, double diag, double* a, int lda,
                        int thread_id, int stream_id);
template void laset_gpu(int m, int n, cuComplex offdiag, cuComplex diag, cuComplex* a, int lda,
                        int thread_id, int stream_id);
template void laset_gpu(int m, int n, cuDoubleComplex offdiag, cuDoubleComplex diag,
                        cuDoubleComplex* a, int lda, int thread_id, int stream_id);

}  // lapack
}  // linalg
}  // dca
