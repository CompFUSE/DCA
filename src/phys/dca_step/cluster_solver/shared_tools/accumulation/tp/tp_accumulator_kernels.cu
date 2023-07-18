// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Weile Wei (wwei9@lsu.ed)
//         Peter Doak (doakpw@ornl.gov)
// Implements the GPU kernels used by the tp_accumulator_gpu

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <complex>

#include "dca/platform/dca_gpu.h"

#include "dca/parallel/util/get_workload.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/type_help.hpp"
#include "dca/util/type_utils.hpp"
#include "dca/linalg/util/atomic_add_cuda.cu.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

using namespace linalg;
using dca::util::ComplexAlias;
using dca::util::castGPUType;
using dca::util::RealAlias;
using phys::FourPointType;
using dca::util::SignType;

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
                                         int nw_freq, const Real beta, int spin) {
  // Computes G = -G0(w1) * M(w1, w2) * G(w2) + (w1 == w2) * beta * G0(w1).

  const int n_rows = nk * nw_freq;
  const int n_cols = n_rows;
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

  const CudaComplex<Real> G0_w1 = G0[k1 + nk * w1];
  const CudaComplex<Real> G0_w2 = G0[k2 + nk * w2];

  const CudaComplex<Real> M_val  = G[id_i + ldg * id_j];
  
  G[id_i + ldg * id_j] = -G0_w1 * M_val * G0_w2;
  if (k1 == k2 && w1 == w2) {
    G[id_i + ldg * id_j] += G0_w1 * beta;
  }

  printf("%f %f %f %f %f %f -- %d %d %d %d %f,%f\n", M_val, G0_w1, G0_w2, spin, k1, k2, w1, w2, G[id_i + ldg * id_j].x, G[id_i + ldg * id_j].y);
}

template <typename Real>
void computeGSingleband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int nk,
                        int nw_freq, const Real beta, cudaStream_t stream, int spin) {
  const int n_rows = nk * nw_freq;
  auto blocks = getBlockSize(n_rows, n_rows);

  computeGSinglebandKernel<<<blocks[0], blocks[1], 0, stream>>>(castGPUType(G), ldg,
                                                                castGPUType(G0), nk, nw_freq, beta, spin);
}

template <typename Real>
__global__ void computeGMultibandKernel(CudaComplex<Real>* __restrict__ G, int ldg,
                                        const CudaComplex<Real>* __restrict__ G0, int ldg0, int nb,
                                        int nk, int nw, Real beta) {
  // Computes G = -G0(w1) * M(w1, w2) * G(w2) + (w1 == w2) * beta * G0(w1).
  // The product is to be intended as matrix-matrix multiplication in band space.

  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;

  if (id_i >= nb * nk * nw || id_j >= nb * nk * nw)
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

  // hmmm... in CPU we now just run over the entire extended range for w1 and w2
  // w1 += nw_pos;

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
    for (int i = 0; i < nb; ++i) {
      const CudaComplex<Real> G_band = G0_w1[b1 + ldg0 * i] * M[i + ldm * j] * G0_w2_val;
      G_val -= G_band;
    }
  }

  if (G0_w1 == G0_w2)
    G_val += G0_w1[b1 + ldg0 * b2] * beta;

  printf("%f %f %f %f %f %f -- %d %d %d %d %d %d %f,%f\n", M[b1 + ldm * b2], G0_w1[b1 + ldg0 * b2], G0_w2[b1 + ldg0 * b2], b1, b2, k1, k2, w1, w2, G_val.x, G_val.y);
}

template <typename Real>
void computeGMultiband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int ldg0,
                       int nb, int nk, int nw, Real beta, cudaStream_t stream) {
  const int n_rows = nb * nk * nw;

  auto get_block_width = [nb] {
    if (nb > 16)
      throw(std::out_of_range("Too many bands."));
    for (int candidate = 16; candidate > 0; --candidate)
      if (!(candidate % nb))
        return candidate;
    return -1;
  };

  const int width = get_block_width();
  // what is the magic 2 for?
  const auto blocks = getBlockSize(n_rows, n_rows, width);

  computeGMultibandKernel<<<blocks[0], blocks[1], width * width * sizeof(std::complex<Real>), stream>>>(
      castGPUType(G), ldg, castGPUType(G0), ldg0, nb, nk, nw, beta);
}

// template <typename Complex>
// __device__ Complex getG(const Complex* __restrict__ G, const int ldg, int k1, int k2, int w1,
//                         int w2, const int b1, const int b2) {
//   const bool is_conj = g4_helper.extendGIndices(k1, k2, w1, w2);

//   const unsigned nb = g4_helper.get_bands();
//   const unsigned nk = g4_helper.get_cluster_size();
//   const unsigned no = nb * nk;

//   unsigned i_idx = b1 + nb * k1 + no * w1;
//   unsigned j_idx = b2 + nb * k2 + no * w2;

//   auto val = G[i_idx + ldg * j_idx];

//   if (!is_conj)
//     return val;
//   else {
//     // For Moire model G_up(-k, -wn) = conj(G_dn(k, wn))
//     // but for a given configuration, SU(2) symmetry is broken by the auxialary spin field
//     // This mean the following code is incorrect and we must extend the calculation to all
//     // frequencies w1, w2, including negative w1.
//     i_idx = 1 - b2 + nb * k1 + no * w1;
//     j_idx = 1 - b1 + nb * k2 + no * w2;
//     val = conj(G[i_idx + ldg * j_idx]);
//     // if (b1==b2)
//     //     return conj(val);
//     // else {
//     //     i_idx = b2 + nb * k1 + no * w1;
//     //     j_idx = b1 + nb * k2 + no * w2;
//     //     val = -conj(G[i_idx + ldg * j_idx]);
//     return val;
//     // }
//   }
//   // return is_conj ? conj(val) : val;
// }

template <typename Scalar, FourPointType type, typename SignType>
__global__ void updateG4Kernel(CudaComplex<RealAlias<Scalar>>* __restrict__ G4,
                               const CudaComplex<RealAlias<Scalar>>* __restrict__ G_up,
                               const int ldgu,
                               const CudaComplex<RealAlias<Scalar>>* __restrict__ G_down,
                               const int ldgd, const SignType factor, const bool atomic,
                               const uint64_t start, const uint64_t end) {
  // TODO: reduce code duplication.
  // TODO: decrease, if possible, register pressure. E.g. a single thread computes all bands.

  const uint64_t local_g4_index =
      static_cast<uint64_t>(blockIdx.x) * static_cast<uint64_t>(blockDim.x) +
      static_cast<uint64_t>(threadIdx.x);

  const uint64_t g4_index = local_g4_index + start;

  if (g4_index >= end) {  // out of domain.
    return;
  }

  Scalar complex_factor;
  dca::linalg::assign(complex_factor, factor);
  const Scalar sign_over_2 = 0.5 * complex_factor;

  int b1, b2, b3, b4, k1, k2, k_ex, w1, w2, w_ex;
  g4_helper.unrollIndex(g4_index, b1, b2, b3, b4, k1, w1, k2, w2, k_ex, w_ex);

  const int nb = g4_helper.get_bands();
  const int nk = g4_helper.get_cluster_size();

  CudaComplex<RealAlias<Scalar>> contribution;
  const unsigned no = nk * nb;
  auto cond_conj = [](const CudaComplex<RealAlias<Scalar>> a, const bool cond) {
    return cond ? conj(a) : a;
  };

  // This code needs to be repeated over and over.  This happens in getGMultiband in the cpu
  // implementation. The gpu code is structed differently so without signficant restructing this
  // can't happen in the extendGIndiciesMultiBand routines.
  auto condSwapAdd = [](int& ia, int& ib, const int ba, const int bb, const bool cond) {
    if (cond) {
      ia += bb;
      ib += ba;
    }
    else {
      ia += ba;
      ib += bb;
    }
  };
  // Compute the contribution to G4. In all the products of Green's function of type Ga * Gb,
  // the dependency on the bands is implied as Ga(b1, b2) * Gb(b2, b3). Sums and differences with
  // the exchange momentum, implies the same operation is performed with the exchange frequency.
  // See tp_accumulator.hpp for more details.
  if constexpr (type == FourPointType::PARTICLE_HOLE_TRANSVERSE) {
    // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, -s)
    int w1_a(w1);
    int w2_a(w2);
    int k1_a(k1);
    int k2_a(k2);
    const bool conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
    int i_a = nb * k1_a + no * w1_a;
    int j_a = nb * k2_a + no * w2_a;
    condSwapAdd(i_a, j_a, b1, b4, conj_a);
    const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
    const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

    int w1_b(g4_helper.addWex(w2, w_ex));
    int w2_b(g4_helper.addWex(w1, w_ex));
    int k1_b = g4_helper.addKex(k2, k_ex);
    int k2_b = g4_helper.addKex(k1, k_ex);
    const bool conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
    int i_b = nb * k1_b + no * w1_b;
    int j_b = nb * k2_b + no * w2_b;
    condSwapAdd(i_b, j_b, b2, b3, conj_b);
    const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
    const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

    contribution = -sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
  }
  else if constexpr (type == FourPointType::PARTICLE_HOLE_MAGNETIC) {
    // The PARTICLE_HOLE_MAGNETIC contribution is computed in two parts:
    // Spin Difference Contribution
    // new scope to reuse local index variables

    // contribution += (\sum_s s * G(k1, k1 + k_ex)) * (\sum_s s * G(k2 + k_ex, k2))
    int k1_a = k1;
    int k2_a = g4_helper.addKex(k1, k_ex);
    int k1_b = g4_helper.addKex(k2, k_ex);
    int k2_b = k2;

    int w1_a(w1);
    int w2_a(g4_helper.addWex(w2, w_ex));
    int w1_b(g4_helper.addWex(w1, w_ex));
    int w2_b(w2);

    // conj_a in this case just tells us whether to swap the band axes additions or not
    bool conj_a = false;
    // if (g4_helper.get_bands() == 1)
    //   conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
    // else
    //   conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
    int i_a = nb * k1_a + no * w1_a;
    int j_a = nb * k2_a + no * w2_a;
    // condSwapAdd(i_a, j_a, b1, b3, conj_a);
    // CudaComplex<RealAlias<Scalar>> Ga =
    //   G_up[i_a + ldgu * j_a] - G_down[i_a + ldgd * j_a];
    // if (i_a == j_a)
    //   Ga += (G_up[i_a + ldgu * j_a] - G_down[i_a + ldgd * j_a]) *

    bool conj_b = false;
    // if (g4_helper.get_bands() == 1)
    //   conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
    // else
    //   conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);
    int i_b = nb * k1_b + no * w1_b;
    int j_b = nb * k2_b + no * w2_b;
    // condSwapAdd(i_b, j_b, b2, b4, conj_b);
    // CudaComplex<RealAlias<Scalar>> Gb =
    //   G_up[i_b + ldgu * j_b] - G_down[i_b + ldgd * j_b];

    // contribution = sign_over_2 * (Ga * Gb);

    // direct contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
    k1_a = k1;
    k2_a = k2;
    k1_b = g4_helper.addKex(k2, k_ex);
    k2_b = g4_helper.addKex(k1, k_ex);

    w1_a = w1;
    w2_a = w2;
    w1_b = g4_helper.addWex(w2, w_ex);
    w2_b = g4_helper.addWex(w1, w_ex);

    if (g4_helper.get_bands() == 1)
      conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
    else
      conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
    i_a = nb * k1_a + no * w1_a;
    j_a = nb * k2_a + no * w2_a;
    i_a += b1;
    j_a += b3;

    CudaComplex<RealAlias<Scalar>> Ga_1 = G_up[i_a + ldgu * j_a];
    CudaComplex<RealAlias<Scalar>> Ga_2 = G_down[i_a + ldgd * j_a];

    if (g4_helper.get_bands() == 1)
      conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
    else
      conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

    i_b = nb * k1_b + no * w1_b;
    j_b = nb * k2_b + no * w2_b;
    i_b += b4;
    j_b += b2;
    CudaComplex<RealAlias<Scalar>> Gb_1 = G_up[i_b + ldgu * j_b];
    CudaComplex<RealAlias<Scalar>> Gb_2 = G_down[i_b + ldgd * j_b];

    contribution = -sign_over_2 * Ga_1 * Gb_1 - sign_over_2 * Ga_2 * Gb_2;
  }
  else if constexpr (type == FourPointType::PARTICLE_HOLE_CHARGE) {
    // The PARTICLE_HOLE_CHARGE contribution is computed in two parts:
    {
      // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);
      bool conj_a = false;
      if (g4_helper.get_bands() == 1)
        conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      else
        conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
      int i_a = nb * k1_a + no * w1_a;
      int j_a = nb * k2_a + no * w2_a;
      if (conj_a) {
        i_a += b4;
        j_a += b1;
      }
      else {
        i_a += b1;
        j_a += b4;
      }

      const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(g4_helper.addWex(w1, w_ex));
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = g4_helper.addKex(k1, k_ex);
      bool conj_b = false;
      if (g4_helper.get_bands() == 1)
        conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      else
        conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

      int i_b = nb * k1_b + no * w1_b;
      int j_b = nb * k2_b + no * w2_b;
      if (conj_b) {
        i_b += b3;
        j_b += b2;
      }
      else {
        i_b += b2;
        j_b += b3;
      }

      const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      contribution = -sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
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
      bool conj_a = false;
      if (g4_helper.get_bands() == 1)
        conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      else
        conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);

      int i_a = nb * k1_a + no * w1_a;
      int j_a = nb * k2_a + no * w2_a;
      if (conj_a) {
        i_a += b1;
        j_a += b2;
      }
      else {
        i_a += b2;
        j_a += b1;
      }

      const CudaComplex<RealAlias<Scalar>> Ga =
          cond_conj(G_up[i_a + ldgu * j_a] + G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(w2);
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = k2;
      bool conj_b = false;
      if (g4_helper.get_bands() == 1)
        conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      else
        conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

      int i_b = nb * k1_b + no * w1_b;
      int j_b = nb * k2_b + no * w2_b;
      if (conj_b) {
        i_b += b4;
        j_b += b3;
      }
      else {
        i_b += b3;
        j_b += b4;
      }

      const CudaComplex<RealAlias<Scalar>> Gb =
          cond_conj(G_up[i_b + ldgu * j_b] + G_down[i_b + ldgd * j_b], conj_b);

      contribution += sign_over_2 * (Ga * Gb);
    }
  }
  else if constexpr (type == FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP) {
    // The PARTICLE_HOLE_LONGITUDINAL_UP_UP contribution is computed in two parts:
    {
      // contribution <- \sum_s G(k1, k1+k_ex, s) * G(k2+k_ex, k2, s)
      int w1_a(w1);
      int w2_a(g4_helper.addWex(w1, w_ex));
      int k1_a = k1;
      int k2_a = g4_helper.addKex(k1, k_ex);
      bool conj_a = false;
      if (g4_helper.get_bands() == 1)
        conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      else
        conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
      int i_a = nb * k1_a + no * w1_a;
      int j_a = nb * k2_a + no * w2_a;
      if (conj_a) {
        i_a += b4;
        j_a += b2;
      }
      else {
        i_a += b2;
        j_a += b4;
      }

      const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(w2);
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = k2;
      bool conj_b = false;
      if (g4_helper.get_bands() == 1)
        conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      else
        conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

      int i_b = nb * k1_b + no * w1_b;
      int j_b = nb * k2_b + no * w2_b;
      if (conj_b) {
        i_b += b1;
        j_b += b3;
      }
      else {
        i_b += b3;
        j_b += b1;
      }

      const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_up[i_b + ldgd * j_b], conj_b);
      const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_down[i_b + ldgu * j_b], conj_b);

      contribution = sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
    }
    {
      // contribution <- -\sum_s G(k1, k2, s) * G(k2 + k_ex, k1 + k_ex, s)
      int w1_a(w1);
      int w2_a(w2);
      int k1_a(k1);
      int k2_a(k2);

      bool conj_a = false;
      if (g4_helper.get_bands() == 1)
        conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
      else
        conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
      int i_a = nb * k1_a + no * w1_a;
      int j_a = nb * k2_a + no * w2_a;
      if (conj_a) {
        i_a += b4;
        j_a += b1;
      }
      else {
        i_a += b1;
        j_a += b4;
      }
      const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

      int w1_b(g4_helper.addWex(w2, w_ex));
      int w2_b(g4_helper.addWex(w1, w_ex));
      int k1_b = g4_helper.addKex(k2, k_ex);
      int k2_b = g4_helper.addKex(k1, k_ex);

      bool conj_b = false;
      if (g4_helper.get_bands() == 1)
        conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
      else
        conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

      int i_b = nb * k1_b + no * w1_b;
      int j_b = nb * k2_b + no * w2_b;
      if (conj_b) {
        i_b += b3;
        j_b += b2;
      }
      else {
        i_b += b2;
        j_b += b3;
      }

      const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_up[i_b + ldgd * j_b], conj_b);
      const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_down[i_b + ldgu * j_b], conj_b);

      contribution += -sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
    }
  }
  else if constexpr (type == FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN) {
    // contribution <- \sum_s G(k1, k1+k_ex, s) * G(k2+k_ex, k2, -s)
    int w1_a(w1);
    int w2_a(g4_helper.addWex(w1, w_ex));
    int k1_a = k1;
    int k2_a = g4_helper.addKex(k1, k_ex);
    bool conj_a = false;
    if (g4_helper.get_bands() == 1)
      conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
    else
      conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
    int i_a = nb * k1_a + no * w1_a;
    int j_a = nb * k2_a + no * w2_a;
    if (conj_a) {
      i_a += b4;
      j_a += b2;
    }
    else {
      i_a += b2;
      j_a += b4;
    }

    const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
    const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

    int w1_b(g4_helper.addWex(w2, w_ex));
    int w2_b(w2);
    int k1_b = g4_helper.addKex(k2, k_ex);
    int k2_b = k2;
    bool conj_b = false;
    if (g4_helper.get_bands() == 1)
      conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
    else
      conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

    int i_b = nb * k1_b + no * w1_b;
    int j_b = nb * k2_b + no * w2_b;
    if (conj_b) {
      i_b += b1;
      j_b += b3;
    }
    else {
      i_b += b3;
      j_b += b1;
    }

    const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
    const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

    contribution = sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
  }
  else if constexpr (type == FourPointType::PARTICLE_PARTICLE_UP_DOWN) {
    // contribution <- -\sum_s G(k_ex - k2, k_ex - k1, s) * G(k2, k1, -s).
    int w1_a(w1);
    int w2_a(w2);
    int k1_a(k1);
    int k2_a(k2);
    bool conj_a = false;
    if (g4_helper.get_bands() == 1)
      conj_a = g4_helper.extendGIndices(k1_a, k2_a, w1_a, w2_a);
    else
      conj_a = g4_helper.extendGIndicesMultiBand(k1_a, k2_a, w1_a, w2_a);
    int i_a = nb * k1_a + no * w1_a;
    int j_a = nb * k2_a + no * w2_a;
    if (conj_a) {
      i_a += b4;
      j_a += b2;
    }
    else {
      i_a += b2;
      j_a += b4;
    }
    const CudaComplex<RealAlias<Scalar>> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
    const CudaComplex<RealAlias<Scalar>> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);

    int w1_b(g4_helper.wexMinus(w1, w_ex));
    int w2_b(g4_helper.wexMinus(w2, w_ex));
    int k1_b = g4_helper.kexMinus(k1, k_ex);
    int k2_b = g4_helper.kexMinus(k2, k_ex);
    bool conj_b = false;
    if (g4_helper.get_bands() == 1)
      conj_b = g4_helper.extendGIndices(k1_b, k2_b, w1_b, w2_b);
    else
      conj_b = g4_helper.extendGIndicesMultiBand(k1_b, k2_b, w1_b, w2_b);

    int i_b = nb * k1_b + no * w1_b;
    int j_b = nb * k2_b + no * w2_b;
    if (conj_b) {
      i_b += b1;
      j_b += b3;
    }
    else {
      i_b += b3;
      j_b += b1;
    }

    const CudaComplex<RealAlias<Scalar>> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);
    const CudaComplex<RealAlias<Scalar>> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

    contribution = sign_over_2 * (Ga_1 * Gb_1 + Ga_2 * Gb_2);
  }

  decltype(G4) const result_ptr = G4 + local_g4_index;
  if (atomic)
    dca::linalg::atomicAdd(result_ptr, contribution);
  else
    *result_ptr += contribution;
}

template <typename Scalar, FourPointType type, typename SignType>
double updateG4(Scalar* G4, const Scalar* G_up, const int ldgu, const Scalar* G_down,
                const int ldgd, const SignType factor, bool atomic, cudaStream_t stream,
                std::size_t start, std::size_t end) {
  constexpr const std::size_t n_threads = 256;
  const unsigned n_blocks = dca::util::ceilDiv(end - start, n_threads);

  using dca::util::GPUTypeConversion;
  updateG4Kernel<dca::util::CUDATypeMap<Scalar>, type><<<n_blocks, n_threads, 0, stream>>>(
      castGPUType(G4), castGPUType(G_up), ldgu, castGPUType(G_down), ldgd,
      GPUTypeConversion(factor), atomic, start, end);

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
    case FourPointType::PARTICLE_HOLE_TRANSVERSE:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    case FourPointType::PARTICLE_HOLE_MAGNETIC:
      // Each update of a G4 entry involves 3 complex additions and 3 complex multiplications.
      return 26. * n_updates;
    case FourPointType::PARTICLE_HOLE_CHARGE:
      // Each update of a G4 entry involves 3 complex additions and 3 complex multiplications.
      return 26. * n_updates;
    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      // Each update of a G4 entry involves 3 complex additions and 4 complex multiplications.
      return 32 * n_updates;
    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    case FourPointType::PARTICLE_PARTICLE_UP_DOWN:
      // Each update of a G4 entry involves 2 complex additions and 2 complex multiplications.
      return 18. * n_updates;
    default:
      throw(std::logic_error("Invalid mode"));
  }
}

// Explicit instantiation.
template void computeGSingleband<float>(std::complex<float>* G, int ldg,
                                        const std::complex<float>* G0, int nk, int nw,
                                        const float beta, cudaStream_t stream, int spin);
template void computeGMultiband<float>(std::complex<float>* G, int ldg,
                                       const std::complex<float>* G0, int ldg0, int nb, int nk,
                                       int nw, float beta, cudaStream_t stream);

template void computeGSingleband<double>(std::complex<double>* G, int ldg,
                                         const std::complex<double>* G0, int nk, int nw_pos,
                                         const double beta, cudaStream_t stream, int spin);
template void computeGMultiband<double>(std::complex<double>* G, int ldg,
                                        const std::complex<double>* G0, int ldg0, int nb, int nk,
                                        int nw_pos, double beta, cudaStream_t stream);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_TRANSVERSE, std::int8_t>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_MAGNETIC, std::int8_t>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_CHARGE, std::int8_t>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                         std::int8_t>(std::complex<float>* G4, const std::complex<float>* G_up,
                                      const int ldgu, const std::complex<float>* G_down,
                                      const int ldgd, const std::int8_t factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
                         std::int8_t>(std::complex<float>* G4, const std::complex<float>* G_up,
                                      const int ldgu, const std::complex<float>* G_down,
                                      const int ldgd, const std::int8_t factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<float>, FourPointType::PARTICLE_PARTICLE_UP_DOWN, std::int8_t>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_TRANSVERSE, std::int8_t>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_MAGNETIC, std::int8_t>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_CHARGE, std::int8_t>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                         std::int8_t>(std::complex<double>* G4, const std::complex<double>* G_up,
                                      const int ldgu, const std::complex<double>* G_down,
                                      const int ldgd, const std::int8_t factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
                         std::int8_t>(std::complex<double>* G4, const std::complex<double>* G_up,
                                      const int ldgu, const std::complex<double>* G_down,
                                      const int ldgd, const std::int8_t factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template double updateG4<std::complex<double>, FourPointType::PARTICLE_PARTICLE_UP_DOWN, std::int8_t>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
    cudaStream_t stream, std::size_t start, std::size_t end);

// complex g0

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_TRANSVERSE, std::complex<float>>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::complex<float> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_MAGNETIC, std::complex<float>>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::complex<float> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_CHARGE, std::complex<float>>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::complex<float> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                std::complex<float>>(std::complex<float>* G4, const std::complex<float>* G_up,
                                     const int ldgu, const std::complex<float>* G_down,
                                     const int ldgd, const std::complex<float> factor, bool atomic,
                                     cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
                std::complex<float>>(std::complex<float>* G4, const std::complex<float>* G_up,
                                     const int ldgu, const std::complex<float>* G_down,
                                     const int ldgd, const std::complex<float> factor, bool atomic,
                                     cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<float>, FourPointType::PARTICLE_PARTICLE_UP_DOWN, std::complex<float>>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const std::complex<float> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_TRANSVERSE, std::complex<double>>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::complex<double> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_MAGNETIC, std::complex<double>>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::complex<double> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_CHARGE, std::complex<double>>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::complex<double> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                std::complex<double>>(std::complex<double>* G4, const std::complex<double>* G_up,
                                      const int ldgu, const std::complex<double>* G_down,
                                      const int ldgd, const std::complex<double> factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
                std::complex<double>>(std::complex<double>* G4, const std::complex<double>* G_up,
                                      const int ldgu, const std::complex<double>* G_down,
                                      const int ldgd, const std::complex<double> factor, bool atomic,
                                      cudaStream_t stream, std::size_t start, std::size_t end);

template <>
double updateG4<std::complex<double>, FourPointType::PARTICLE_PARTICLE_UP_DOWN, std::complex<double>>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const std::complex<double> factor,
    bool atomic, cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_TRANSVERSE>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_MAGNETIC>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_CHARGE>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_PARTICLE_UP_DOWN>(
//   std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
//   const std::complex<float>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_TRANSVERSE>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_MAGNETIC>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_CHARGE>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

// template<> double updateG4< FourPointType::PARTICLE_PARTICLE_UP_DOWN>(
//   std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
//   const std::complex<double>* G_down, const int ldgd, const std::int8_t factor, bool atomic,
//   cudaStream_t stream, std::size_t start, std::size_t end);

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca
