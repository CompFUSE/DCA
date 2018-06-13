// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/util/cuda_blocks.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {

__global__ void setRightSectorToIdKernel(double* m, const int ldm, const int n0, const int n_max) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y + n0;

  if (i >= n_max || j >= n_max)
    return;

  m[i + ldm * j] = (i == j) ? 1. : 0.;
}

void setRightSectorToId(double* m, const int ldm, const int n0, const int n_max, cudaStream_t stream) {
  auto blocks = dca::util::getBlockSize(n_max, n_max - n0);

  setRightSectorToIdKernel<<<blocks[0], blocks[1], 0, stream>>>(m, ldm, n0, n_max);
}

__global__ void computeGLeftKernel(MatrixView G, const MatrixView M, const double* f, int n_init) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= G.nrRows() || j >= n_init)
    return;

  G(i, j) = (M(i, j) * f[j] - double(i == j)) / (f[j] - 1);
}

void computeGLeft(MatrixView& G, const MatrixView& M, const double* f, int n_init,
                  cudaStream_t stream) {
  if (n_init == 0)
    return;
  const int n = G.nrRows();
  const auto blocks = dca::util::getBlockSize(n, n_init);

  computeGLeftKernel<<<blocks[0], blocks[1], 0, stream>>>(G, M, f, n_init);
}

__global__ void multiplyByFFactorKernel(MatrixView M, const double* f_vals,
                                        const bool inverse_factor, const bool row_factor) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= M.nrRows() || j >= M.nrCols())
    return;

  double factor;
  if (row_factor)
    factor = -(f_vals[i] - 1.);
  else
    factor = f_vals[j] - 1.;

  if (inverse_factor)
    M(i, j) /= factor;
  else
    M(i, j) *= factor;
}

void multiplyByFFactor(MatrixView& M, const double* f_vals, bool inverse_factor, bool row_factor,
                       cudaStream_t stream) {
  if (M.nrCols() == 0 || M.nrRows() == 0)
    return;
  const auto blocks = dca::util::getBlockSize(M.nrRows(), M.nrCols());

  multiplyByFFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(M, f_vals, inverse_factor, row_factor);
}

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca
