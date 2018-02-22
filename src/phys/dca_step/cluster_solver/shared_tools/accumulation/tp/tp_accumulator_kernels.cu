// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implements the GPU kernels used by the DFT algorithm.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
namespace global {
// dca::phys::solver::accumulator::details::global::
G4HelperManager helper;
}  // global
// dca::phys::solver::accumulator::details::

using namespace linalg;
using linalg::util::CudaComplex;
using linalg::util::castCudaComplex;

std::array<dim3, 2> getBlockSize(const uint i, const uint j, const uint block_size = 32) {
  const uint n_threads_i = std::min(block_size, i);
  const uint n_threads_j = std::min(block_size, j);
  if (n_threads_i * n_threads_j > 32 * 32)
    throw(std::logic_error("Block size is too big"));

  const uint n_blocks_i = dca::util::ceilDiv(i, n_threads_i);
  const uint n_blocks_j = dca::util::ceilDiv(j, n_threads_j);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j), dim3(n_threads_i, n_threads_j)};
}

template <typename Real>
__global__ void computeGSinglebandKernel(CudaComplex<Real>* __restrict__ G, int ldg,
                                         const CudaComplex<Real>* __restrict__ G0, int nk,
                                         int nw_pos, const Real beta) {
  const int n_rows = nk * nw_pos;
  const int n_cols = n_rows * 2;
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= n_rows || id_j >= n_cols)
    return;

  auto get_indices = [=](const int id, int& k, int& w) {
    w = id / nk;
    k = id - nk * w;
  };
  int w1, w2, k1, k2;
  get_indices(id_i, k1, w1);
  get_indices(id_j, k2, w2);

  const CudaComplex<Real> G0_w1 = G0[k1 + nk * (w1 + nw_pos)];
  const CudaComplex<Real> G0_w2 = G0[k2 + nk * w2];

  G[id_i + ldg * id_j] *= -G0_w1 * G0_w2;
  if (k1 == k2 && w1 + nw_pos == w2) {
    G[id_i + ldg * id_j] += G0_w1 * beta;
  }
}

template <typename Real>
void computeGSingleband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int nk,
                        int nw_pos, const Real beta, cudaStream_t stream) {
  const int n_rows = nk * nw_pos;
  auto blocks = getBlockSize(n_rows, n_rows * 2);

  computeGSinglebandKernel<<<blocks[0], blocks[1], 0, stream>>>(
      castCudaComplex(G), ldg, castCudaComplex(G0), nk, nw_pos, beta);
}

template <typename Real>
__global__ void computeGMultibandKernel(CudaComplex<Real>* __restrict__ G, int ldg,
                                        const CudaComplex<Real>* __restrict__ G0, int ldg0, int nb,
                                        int nk, int nw_pos, Real beta) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;

  assert(id_i < nb * nk * nw_pos);
  assert(id_j < nb * nk * nw_pos * 2);

  const int no = nb * nk;
  auto get_indices = [=](int id, int& b, int& k, int& w) {
    w = id / no;
    id -= w * no;
    k = id / nb;
    b = id - k * nb;
  };
  int w1, w2, k1, k2, b1, b2;
  get_indices(id_i, b1, k1, w1);
  get_indices(id_j, b2, k2, w2);
  w1 += nw_pos;

  // Note: cuda does not support templated shared memory.
  extern __shared__ char shared_mem[];
  CudaComplex<Real>* const M_block = reinterpret_cast<CudaComplex<Real>*>(shared_mem);
  const int local_row_start = (threadIdx.y / nb) * nb;
  const int local_col_start = (threadIdx.x / nb) * nb;
  const int ldm = blockDim.y;
  CudaComplex<Real>* const M = M_block + local_row_start + ldm * local_col_start;

  CudaComplex<Real>& G_val = G[id_i + ldg * id_j];
  M[b1 + ldm * b2] = G_val;
  __syncthreads();

  const CudaComplex<Real>* const G0_w1 = G0 + nb * k1 + no * w1;
  const CudaComplex<Real>* const G0_w2 = G0 + nb * k2 + no * w2;

  G_val.x = G_val.y = 0;
  for (int j = 0; j < nb; ++j) {
    const CudaComplex<Real> G0_w2_val = G0_w2[j + ldg0 * b2];
    for (int i = 0; i < nb; ++i)
      G_val -= G0_w1[b1 + ldg0 * i] * M[i + ldm * j] * G0_w2_val;
  }

  if (G0_w1 == G0_w2)
    G_val += G0_w1[b1 + ldg0 * b2] * beta;
}

template <typename Real>
void computeGMultiband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int ldg0,
                       int nb, int nk, int nw_pos, Real beta, cudaStream_t stream) {
  const int n_rows = nb * nk * nw_pos;

  auto get_block_width = [nb] {
    if (nb > 16)
      throw(std::logic_error("Too many bands."));
    for (int candidate = 16; candidate > 0; --candidate)
      if (!(candidate % nb))
        return candidate;
    return -1;
  };
  const static int width = get_block_width();

  const static auto blocks = getBlockSize(n_rows, n_rows * 2, width);

  computeGMultibandKernel<<<blocks[0], blocks[1], width * width * sizeof(std::complex<Real>), stream>>>(
      castCudaComplex(G), ldg, castCudaComplex(G0), ldg0, nb, nk, nw_pos, beta);
}

void initializeG4Helpers(int nb, int nk, int nw_pos, int delta_k, int delta_w, const int* add_k,
                         int lda, const int* sub_k, int lds) {
  if (!global::helper.isInitialized())
    global::helper.set(nb, nk, nw_pos, delta_k, delta_w, add_k, lda, sub_k, lds);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

// Include specializations for each mode
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/modes/particle_particle_up_down.inc"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/modes/particle_hole_transverse.inc"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/modes/particle_hole_charge.inc"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/modes/particle_hole_magnetic.inc"

// Explicit instantiation.
template void computeGSingleband<float>(std::complex<float>* G, int ldg,
                                        const std::complex<float>* G0, int nk, int nw,
                                        const float beta, cudaStream_t stream);
template void computeGMultiband<float>(std::complex<float>* G, int ldg,
                                       const std::complex<float>* G0, int ldg0, int nb, int nk,
                                       int nw, float beta, cudaStream_t stream);

template void computeGSingleband<double>(std::complex<double>* G, int ldg,
                                         const std::complex<double>* G0, int nk, int nw_pos,
                                         const double beta, cudaStream_t stream);
template void computeGMultiband<double>(std::complex<double>* G, int ldg,
                                        const std::complex<double>* G0, int ldg0, int nb, int nk,
                                        int nw_pos, double beta, cudaStream_t stream);

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca
