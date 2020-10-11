// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of the device methods used by the CT-INT submatrix walker.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <complex>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/util/cuda_blocks.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {

using namespace dca::linalg;

template <typename Scalar>
__global__ void setRightSectorToIdKernel(Scalar* m, const int ldm, const int n0, const int n_max) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y + n0;

  if (i >= n_max || j >= n_max)
    return;

  using Real = dca::util::Real<Scalar>;
  assign(m[i + ldm * j], (i == j) ? Real(1.) : Real(0.));
}

template <typename Scalar>
void setRightSectorToId(Scalar* m, const int ldm, const int n0, const int n_max, cudaStream_t stream) {
  auto blocks = dca::util::getBlockSize(n_max, n_max - n0);

  setRightSectorToIdKernel<<<blocks[0], blocks[1], 0, stream>>>(castCuda(m), ldm, n0, n_max);
}

template void setRightSectorToId(float*, const int, const int, const int, cudaStream_t);
template void setRightSectorToId(double*, const int, const int, const int, cudaStream_t);
template void setRightSectorToId(std::complex<double>*, const int, const int, const int,
                                 cudaStream_t);
template void setRightSectorToId(std::complex<float>*, const int, const int, const int, cudaStream_t);

template <typename Scalar, typename Real>
__global__ void computeGLeftKernel(MatrixView<Scalar, GPU> G, const MatrixView<Scalar, GPU> M,
                                   const Real* __restrict__ f, int n_init) {
  const int i_t = threadIdx.x + blockDim.x * blockIdx.x;
  const int stride = blockDim.x * gridDim.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (j >= n_init)
    return;

  const Real factor = 1. / (f[j] - 1);
  const Real fj = f[j];

  for (int i = i_t; i < G.nrRows(); i += stride)
    G(i, j) = (M(i, j) * fj - Real(i == j)) * factor;
}

template <typename Scalar, typename Real>
void computeGLeft(MatrixView<Scalar, GPU>& G, const MatrixView<Scalar, GPU>& M, const Real* f,
                  int n_init, cudaStream_t stream) {
  if (n_init == 0)
    return;
  const int n = G.nrRows();

  constexpr int thread_j = 4;
  constexpr int thread_i = 64;
  dim3 threads(thread_i, thread_j);
  dim3 blocks(std::max(n / (10 * thread_i), 1), dca::util::ceilDiv(n_init, thread_j));

  computeGLeftKernel<<<blocks, threads, 0, stream>>>(castCuda(G), castCuda(M), f, n_init);
}

template void computeGLeft(MatrixView<float, GPU>&, const MatrixView<float, GPU>&, const float*,
                           int, cudaStream_t);
template void computeGLeft(MatrixView<std::complex<float>, GPU>&,
                           const MatrixView<std::complex<float>, GPU>&, const float*, int,
                           cudaStream_t);
template void computeGLeft(MatrixView<double, GPU>&, const MatrixView<double, GPU>&, const double*,
                           int, cudaStream_t);
template void computeGLeft<std::complex<double>>(MatrixView<std::complex<double>, GPU>&,
                                                 const MatrixView<std::complex<double>, GPU>&,
                                                 const double*, int, cudaStream_t);

template <typename Scalar, typename Real>
__global__ void multiplyByFColFactorKernel(MatrixView<Scalar, GPU> M, const Real* f_vals) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= M.nrRows() || j >= M.nrCols())
    return;

  const Real factor = f_vals[j] - 1.;
  M(i, j) *= factor;
}

template <typename Scalar, typename Real>
void multiplyByFColFactor(MatrixView<Scalar, GPU>& M, const Real* f_vals, cudaStream_t stream) {
  if (M.nrCols() == 0 || M.nrRows() == 0)
    return;
  const auto blocks = dca::util::getBlockSize(M.nrRows(), M.nrCols());

  multiplyByFColFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(castCuda(M), f_vals);
}

template void multiplyByFColFactor(MatrixView<float, GPU>&, const float*, cudaStream_t);
template void multiplyByFColFactor(MatrixView<std::complex<float>, GPU>&, const float*, cudaStream_t);
template void multiplyByFColFactor(MatrixView<double, GPU>&, const double*, cudaStream_t);
template void multiplyByFColFactor<std::complex<double>>(MatrixView<std::complex<double>, GPU>&,
                                                         const double*, cudaStream_t);

template <typename Scalar, typename Real>
__global__ void multiplyByInverseFFactorKernel(const MatrixView<Scalar, GPU> m_in,
                                               MatrixView<Scalar, GPU> m_out, const Real* f_vals) {
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= m_in.nrRows() || j >= m_in.nrCols())
    return;

  const Real factor = -(f_vals[i] - 1.);
  m_out(i, j) = m_in(i, j) * factor;
}

template <typename Scalar, typename Real>
void multiplyByInverseFFactor(const MatrixView<Scalar, GPU>& m_in, MatrixView<Scalar, GPU>& m_out,
                              const Real* f_vals, cudaStream_t stream) {
  assert(m_in.nrRows() == m_out.nrRows() && m_in.nrCols() == m_out.nrCols());
  if (m_in.nrCols() == 0 || m_in.nrRows() == 0)
    return;
  const auto blocks = dca::util::getBlockSize(m_in.nrRows(), m_out.nrCols());

  using dca::linalg::castCuda;
  multiplyByInverseFFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(castCuda(m_in),
                                                                      castCuda(m_out), f_vals);
}

template void multiplyByInverseFFactor(const MatrixView<float, GPU>&, MatrixView<float, GPU>&,
                                       const float*, cudaStream_t);
template void multiplyByInverseFFactor(const MatrixView<std::complex<float>, GPU>&,
                                       MatrixView<std::complex<float>, GPU>&, const float*,
                                       cudaStream_t);
template void multiplyByInverseFFactor(const MatrixView<double, GPU>&, MatrixView<double, GPU>&,
                                       const double*, cudaStream_t);
template void multiplyByInverseFFactor(const MatrixView<std::complex<double>, GPU>&,
                                       MatrixView<std::complex<double>, GPU>&, const double*,
                                       cudaStream_t);

template <typename Scalar, typename Real>
__global__ void divideByGammaFactorKernel(MatrixView<Scalar, GPU> m,
                                          const std::pair<int, Real>* gamma_indices,
                                          const int n_indices) {
  // TODO: loop over a number of j indices.
  const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int j = threadIdx.y + blockDim.y * blockIdx.y;
  if (i >= n_indices || j >= m.nrCols())
    return;

  const int p = gamma_indices[i].first;
  assert(p < m.nrRows());

  m(p, j) /= Real(1.) + gamma_indices[i].second;
}

template <typename Scalar, typename Real>
void divideByGammaFactor(MatrixView<Scalar, GPU> m, const std::pair<int, Real>* gamma_indices,
                         const int n_indices, cudaStream_t stream) {
  const auto blocks = dca::util::getBlockSize(n_indices, m.nrCols());

  divideByGammaFactorKernel<<<blocks[0], blocks[1], 0, stream>>>(castCuda(m), gamma_indices,
                                                                 n_indices);
}

template void divideByGammaFactor(MatrixView<float, GPU>, const std::pair<int, float>*, const int,
                                  cudaStream_t);
template void divideByGammaFactor(MatrixView<std::complex<float>, GPU>,
                                  const std::pair<int, float>*, const int, cudaStream_t);
template void divideByGammaFactor(MatrixView<double, GPU>, const std::pair<int, double>*, const int,
                                  cudaStream_t);
template void divideByGammaFactor(MatrixView<std::complex<double>, GPU>,
                                  const std::pair<int, double>*, const int, cudaStream_t);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
