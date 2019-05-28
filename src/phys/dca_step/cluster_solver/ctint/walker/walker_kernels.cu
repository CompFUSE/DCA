// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

__global__ void computeGLeftKernel(MatrixView G, const MatrixView M, const double* __restrict__ f, int n_init) {
  const int i_t = threadIdx.x + blockDim.x * blockIdx.x;
  const int stride = blockDim.x * gridDim.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (j >= n_init)
    return;

  const double factor = 1. / (f[j] - 1);
  const double fj = f[j];

  for (int i = i_t; i < G.nrRows(); i += stride)
    G(i, j) = (M(i, j) * fj - double(i == j)) * factor;
}

void computeGLeft(MatrixView& G, const MatrixView& M, const double* f, int n_init,
                  cudaStream_t stream) {
  if (n_init == 0)
    return;
  const int n = G.nrRows();

  constexpr int thread_j = 4;
  constexpr int thread_i = 64;
  dim3 threads(thread_i, thread_j);
  dim3 blocks(std::max(n / (10 * thread_i), 1), util::ceilDiv(n_init, thread_j));

  computeGLeftKernel<<<blocks, threads, 0, stream>>>(G, M, f, n_init);
}

__global__ void multiplyByFColFactorKernel(MatrixView M, const double* f_vals) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= M.nrRows() || j >= M.nrCols())
    return;

  const double factor = f_vals[j] - 1.;
  M(i, j) *= factor;
}

void multiplyByFColFactor(MatrixView& M, const double* f_vals, cudaStream_t stream) {
  if (M.nrCols() == 0 || M.nrRows() == 0)
    return;
  const auto blocks = dca::util::getBlockSize(M.nrRows(), M.nrCols());

  multiplyByFColFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(M, f_vals);
}

__global__ void multiplyByInverseFFactorKernel(const MatrixView m_in, MatrixView m_out,
                                               const double* f_vals) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= m_in.nrRows() || j >= m_in.nrCols())
    return;

  const double factor = -(f_vals[i] - 1.);
  m_out(i, j) = factor * m_in(i, j);
}

void multiplyByInverseFFactor(const MatrixView& m_in, MatrixView& m_out, const double* f_vals,
                              cudaStream_t stream) {
  assert(m_in.nrRows() == m_out.nrRows() && m_in.nrCols() == m_out.nrCols());
  if (m_in.nrCols() == 0 || m_in.nrRows() == 0)
    return;
  const auto blocks = dca::util::getBlockSize(m_in.nrRows(), m_out.nrCols());

  multiplyByInverseFFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(m_in, m_out, f_vals);
}

__global__ void divideByGammaFactorKernel(MatrixView m, const std::pair<int, double>* gamma_indices,
                                          const int n_indices) {
  // TODO: loop over a number of j indices.
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= n_indices || j >= m.nrCols())
    return;

  const int p = gamma_indices[i].first;
  assert(p < m.nrRows());

  m(p, j) /= 1. + gamma_indices[i].second;
}

void divideByGammaFactor(MatrixView m, const std::pair<int, double>* gamma_indices,
                         const int n_indices, cudaStream_t stream) {
  const auto blocks = dca::util::getBlockSize(n_indices, m.nrCols());

  divideByGammaFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(m, gamma_indices, n_indices);
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
