// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Weile Wei (wwei9@lsu.ed)
//
// Implements the GPU kernels used by the DFT algorithm.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <complex>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/parallel/util/get_workload.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/atomic_add_cuda.cu.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
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
  // Computes G = -G0(w1) * M(w1, w2) * G(w2) + (w1 == w2) * beta * G0(w1).

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
  // Computes G = -G0(w1) * M(w1, w2) * G(w2) + (w1 == w2) * beta * G0(w1).
  // The product is to be intended as matrix-matrix multiplication in band space.

  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;

  if (id_i >= nb * nk * nw_pos || id_j >= nb * nk * nw_pos * 2)
    return;

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
      throw(std::out_of_range("Too many bands."));
    for (int candidate = 16; candidate > 0; --candidate)
      if (!(candidate % nb))
        return candidate;
    return -1;
  };

  const int width = get_block_width();
  const auto blocks = getBlockSize(n_rows, n_rows * 2, width);

  computeGMultibandKernel<<<blocks[0], blocks[1], width * width * sizeof(std::complex<Real>), stream>>>(
      castCudaComplex(G), ldg, castCudaComplex(G0), ldg0, nb, nk, nw_pos, beta);
}

template <typename Real, FourPointType type>
__global__ void updateG4Kernel(CudaComplex<Real>* __restrict__ G4,
                               const CudaComplex<Real>* __restrict__ G_up, const int ldgu,
                               const CudaComplex<Real>* __restrict__ G_down, const int ldgd,
                               const int sign, const bool atomic, const uint64_t start,
                               const uint64_t end) {
  // TODO: reduce code duplication.
  // TODO: decrease, if possible, register pressure. E.g. a single thread computes all bands.

  const uint64_t local_g4_index =
      static_cast<uint64_t>(blockIdx.x) * static_cast<uint64_t>(blockDim.x) +
      static_cast<uint64_t>(threadIdx.x);

  const uint64_t g4_index = local_g4_index + start;

  if (g4_index >= end) {  // out of domain.
    return;
  }

  unsigned b1, b2, b3, b4, k1, k2, k_ex, w1, w2, w_ex;
  g4_helper.unrollIndex(g4_index, b1, b2, b3, b4, k1, w1, k2, w2, k_ex, w_ex);

  const unsigned nb = g4_helper.get_bands();
  const unsigned nk = g4_helper.get_cluster_size();

  CudaComplex<Real> contribution;
  const unsigned no = nk * nb;
  auto cond_conj = [](const CudaComplex<Real> a, const bool cond) { return cond ? conj(a) : a; };

  // Compute the contribution to G4. In all the products of Green's function of type Ga * Gb,
  // the dependency on the bands is implied as Ga(b1, b2) * Gb(b2, b3). Sums and differences with
  // the exchange momentum, implies the same operation is performed with the exchange frequency.
  // See tp_accumulator.hpp for more details.
  switch (type) {
    case PARTICLE_HOLE_TRANSVERSE: {
      // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, -s)
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b4 + nb * k2_a + no * w2_a;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(g4_helper.addWex(w1, w_ex));
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = g4_helper.addKex(k1, k_ex);
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b3 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
      const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      contribution = -(Ga_1 * Gb_1 + Ga_2 * Gb_2);
    } break;

    // The PARTICLE_HOLE_MAGNETIC contribution is computed in two parts:
    case PARTICLE_HOLE_MAGNETIC: {
      // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b4 + nb * k2_a + no * w2_a;
      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(g4_helper.addWex(w1, w_ex));
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = g4_helper.addKex(k1, k_ex);
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b3 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      contribution = -(Ga_1 * Gb_1 + Ga_2 * Gb_2);
    }
      // Spin Difference Contribution
      // new scope to reuse local index variables
      {
        // contribution += (\sum_s s * G(k1, k1 + k_ex)) * (\sum_s s * G(k2 + k_ex, k2))
        int w1_a(w1);
        int w2_a(g4_helper.addWex(w1, w_ex));
        int k1_a = k1;
        int k2_a = g4_helper.addKex(k1, k_ex);
        const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
        const int i_a = b1 + nb * k1_a + no * w1_a;
        const int j_a = b3 + nb * k2_a + no * w2_a;

        const CudaComplex<Real> Ga =
            cond_conj(G_up[i_a + ldgu * j_a] - G_down[i_a + ldgd * j_a], conj_a);

        int w1_b(g4_helper.addWex(w2, w_ex));
        int w2_b(w2);
        int k1_b = g4_helper.addKex(k2, k_ex);
        int k2_b = k2;
        const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
        const int i_b = b2 + nb * k1_b + no * w1_b;
        const int j_b = b4 + nb * k2_b + no * w2_b;

        const CudaComplex<Real> Gb =
            cond_conj(G_up[i_b + ldgu * j_b] - G_down[i_b + ldgd * j_b], conj_b);

        contribution += (Ga * Gb);
      }
      break;

    // The PARTICLE_HOLE_CHARGE contribution is computed in two parts:
    case PARTICLE_HOLE_CHARGE: {
      // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b4 + nb * k2_a + no * w2_a;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(g4_helper.addWex(w1, w_ex));
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = g4_helper.addKex(k1, k_ex);
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b3 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      contribution = -(Ga_1 * Gb_1 + Ga_2 * Gb_2);
    }
      // Spin Difference Contribution
      // new scope to reuse local index variables
      {
        // contribution += (\sum_s G(k1, k1 + k_ex, s)) * (\sum_s G(k2 + k_ex, k2, s))
        // TODO: pull into function, index setting code is identical for Spin cases
        int w1_a(w1);
        int w2_a(g4_helper.addWex(w1, w_ex));
        int k1_a = k1;
        int k2_a = g4_helper.addKex(k1, k_ex);
        const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
        const int i_a = b1 + nb * k1_a + no * w1_a;
        const int j_a = b3 + nb * k2_a + no * w2_a;

        const CudaComplex<Real> Ga =
            cond_conj(G_up[i_a + ldgu * j_a] + G_down[i_a + ldgd * j_a], conj_a);

        int w1_b(g4_helper.addWex(w2, w_ex));
        int w2_b(w2);
        int k1_b = g4_helper.addKex(k2, k_ex);
        int k2_b = k2;
        const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
        const int i_b = b2 + nb * k1_b + no * w1_b;
        const int j_b = b4 + nb * k2_b + no * w2_b;

        const CudaComplex<Real> Gb =
            cond_conj(G_up[i_b + ldgu * j_b] + G_down[i_b + ldgd * j_b], conj_b);

        contribution += (Ga * Gb);
      }
      break;

      // The PARTICLE_HOLE_LONGITUDINAL_UP_UP contribution is computed in two parts:
    case PARTICLE_HOLE_LONGITUDINAL_UP_UP: {
      // contribution <- \sum_s G(k1, k1+k_ex, s) * G(k2+k_ex, k2, s)
      int w1_a(w1);
      int w2_a(g4_helper.addWex(w1, w_ex));
      int k1_a = k1;
      int k2_a = g4_helper.addKex(k1, k_ex);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b3 + nb * k2_a + no * w2_a;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(w2);
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = k2;
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b4 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgd * j_b], conj_b);
      const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgu * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2);
    }
      {
        // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
        int w1_a(w1);
        int w2_a(w2);
        int k1_a(k1);
        int k2_a(k2);
        const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
        const int i_a = b1 + nb * k1_a + no * w1_a;
        const int j_a = b4 + nb * k2_a + no * w2_a;

        const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
        const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

        int w1_b(g4_helper.addWex(w2, w_ex));
        int w2_b(g4_helper.addWex(w1, w_ex));
        int k1_b = g4_helper.addKex(k2, k_ex);
        int k2_b = g4_helper.addKex(k1, k_ex);
        const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
        const int i_b = b2 + nb * k1_b + no * w1_b;
        const int j_b = b3 + nb * k2_b + no * w2_b;

        const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgd * j_b], conj_b);
        const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgu * j_b], conj_b);

        contribution += -(Ga_1 * Gb_1 + Ga_2 * Gb_2);
      }
      break;

    case PARTICLE_HOLE_LONGITUDINAL_UP_DOWN: {
      // contribution <- \sum_s G(k1, k1+k_ex, s) * G(k2+k_ex, k2, -s)
      int w1_a(w1);
      int w2_a(g4_helper.addWex(w1, w_ex));
      int k1_a = k1;
      int k2_a = g4_helper.addKex(k1, k_ex);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b3 + nb * k2_a + no * w2_a;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(w2);
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = k2;
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b4 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
      const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2);
    } break;

    case PARTICLE_PARTICLE_UP_DOWN: {
      // contribution <- -\sum_s G(k_ex - k2, k_ex - k1, s) * G(k2, k1, -s).
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);
      const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      const int i_a = b1 + nb * k1_a + no * w1_a;
      const int j_a = b3 + nb * k2_a + no * w2_a;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.wexMinus(w1, w_ex));
      int w2_b(g4_helper.wexMinus(w2, w_ex));
      int k1_b = g4_helper.kexMinus(k1, k_ex);
      int k2_b = g4_helper.kexMinus(k2, k_ex);
      const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      const int i_b = b2 + nb * k1_b + no * w1_b;
      const int j_b = b4 + nb * k2_b + no * w2_b;

      const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
      const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2);
    } break;
    default:  // abort
      asm("trap;");
  }

  CudaComplex<Real>* const result_ptr = G4 + local_g4_index;
  if (atomic)
    dca::linalg::atomicAdd(result_ptr, contribution * 0.5 * sign);
  else
    *result_ptr += contribution * 0.5 * sign;
}

template <typename Real, FourPointType type>
float updateG4(std::complex<Real>* G4, const std::complex<Real>* G_up, const int ldgu,
               const std::complex<Real>* G_down, const int ldgd, const int sign, bool atomic,
               cudaStream_t stream, std::size_t start, std::size_t end) {
  constexpr const std::size_t n_threads = 256;
  const unsigned n_blocks = dca::util::ceilDiv(end - start, n_threads);

  updateG4Kernel<Real, type>
      <<<n_blocks, n_threads, 0, stream>>>(castCudaComplex(G4), castCudaComplex(G_up), ldgu,
                                           castCudaComplex(G_down), ldgd, sign, atomic, start, end);

  // Check for errors.
  auto err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    linalg::util::printErrorMessage(err, __FUNCTION__, __FILE__, __LINE__);
    throw(std::runtime_error("CUDA failed to launch the G4 kernel."));
  }

  const std::size_t n_updates = end - start;
  switch (type) {
      // Note: sign flips  are ignored and a single complex * real multiplication is
      // present in all modes.
    case PARTICLE_HOLE_TRANSVERSE:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    case PARTICLE_HOLE_MAGNETIC:
      // Each update of a G4 entry involves 3 complex additions and 3 complex multiplications.
      return 26. * n_updates;
    case PARTICLE_HOLE_CHARGE:
      // Each update of a G4 entry involves 3 complex additions and 3 complex multiplications.
      return 26. * n_updates;
    case PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      // Each update of a G4 entry involves 3 complex additions and 4 complex multiplications.
      return 32 * n_updates;
    case PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    case PARTICLE_PARTICLE_UP_DOWN:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    default:
      throw(std::logic_error("Invalid mode"));
  }
}

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

template float updateG4<float, PARTICLE_HOLE_TRANSVERSE>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<float, PARTICLE_HOLE_MAGNETIC>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<float, PARTICLE_HOLE_CHARGE>(std::complex<float>* G4,
                                                     const std::complex<float>* G_up, const int ldgu,
                                                     const std::complex<float>* G_down,
                                                     const int ldgd, const int sign, bool atomic,
                                                     cudaStream_t stream, std::size_t start,
                                                     std::size_t end);

template float updateG4<float, PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<float, PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<float, PARTICLE_PARTICLE_UP_DOWN>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_HOLE_TRANSVERSE>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_HOLE_MAGNETIC>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_HOLE_CHARGE>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template float updateG4<double, PARTICLE_PARTICLE_UP_DOWN>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int sign, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
