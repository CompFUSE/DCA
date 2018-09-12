// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
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
#include "dca/linalg/util/atomic_add_cuda.cu.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/error_cuda.hpp"
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

std::array<dim3, 2> getBlockSize3D(const uint i, const uint j, const uint k) {
  const uint n_threads_k = std::min(uint(8), k);
  const uint max_block_size_ij = n_threads_k > 1 ? 8 : 32;
  const uint n_threads_i = std::min(max_block_size_ij, i);
  const uint n_threads_j = std::min(max_block_size_ij, j);

  const uint n_blocks_i = dca::util::ceilDiv(i, n_threads_i);
  const uint n_blocks_j = dca::util::ceilDiv(j, n_threads_j);
  const uint n_blocks_k = dca::util::ceilDiv(k, n_threads_k);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j, n_blocks_k),
                             dim3(n_threads_i, n_threads_j, n_blocks_k)};
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
  const static int width = get_block_width();

  const static auto blocks = getBlockSize(n_rows, n_rows * 2, width);

  computeGMultibandKernel<<<blocks[0], blocks[1], width * width * sizeof(std::complex<Real>), stream>>>(
      castCudaComplex(G), ldg, castCudaComplex(G0), ldg0, nb, nk, nw_pos, beta);
}

template <typename Real, FourPointType type>
__global__ void updateG4Kernel(CudaComplex<Real>* __restrict__ G4,
                               const CudaComplex<Real>* __restrict__ G_up, const int ldgu,
                               const CudaComplex<Real>* __restrict__ G_down, const int ldgd,
                               const int nb, const int nk, const int nw, const int nw_exchange,
                               const int nk_exchange, const int sign, const G4Helper helper) {
  const int size = nk * nw * nb * nb;
  // id_i is a linearized index of b1, b2, k1, k2.
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  // id_j is a linearized index of b3, b4, k2, k_ex.
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  // id_z is a linearized index of k_ex, w_ex.
  const int id_z = blockIdx.z * blockDim.z + threadIdx.z;
  if (id_i >= size || id_j >= size || id_z >= nw_exchange * nk_exchange)
    return;

  // Unroll id_i and id_j.
  const int step2 = nb * nb;
  const int step1 = step2 * nk;
  auto get_indices = [=](int id, int& b1, int& b2, int& k, int& w) {
    w = id / step1;
    id -= w * step1;
    k = id / step2;
    id -= k * step2;
    b2 = id / nb;
    b1 = id - nb * b2;
  };
  int w1, w2, k1, k2, b1, b2, b3, b4;
  get_indices(id_i, b1, b2, k1, w1);
  get_indices(id_j, b3, b4, k2, w2);

  // Unroll the exchange index id_z = k_ex + nk_exchange * w_ex.
  const int w_ex = id_z / nk_exchange;
  const int k_ex = id_z - w_ex * nk_exchange;

  CudaComplex<Real> contribution;
  const Real factor = 0.5 * sign;
  const int no = nk * nb;
  auto cond_conj = [](const CudaComplex<Real> a, const bool cond) { return cond ? conj(a) : a; };

  // Compute the contribution.
  switch (type) {
    case PARTICLE_HOLE_TRANSVERSE: {
      int w1_a(w1), w2_a(w2);
      const bool conj_a = helper.extendWIndices(w1_a, w2_a);
      const int i_a = b1 + nb * k1 + no * w1_a;
      const int j_a = b4 + nb * k2 + no * w2_a;

      int w1_b(helper.addWex(w2, w_ex));
      int w2_b(helper.addWex(w1, w_ex));
      const bool conj_b = helper.extendWIndices(w1_b, w2_b);
      const int i_b = b2 + nb * helper.addKex(k2, k_ex) + no * w1_b;
      const int j_b = b3 + nb * helper.addKex(k1, k_ex) + no * w2_b;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
      const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2) * (-factor);
    } break;

    case PARTICLE_HOLE_MAGNETIC: {
      int w1_a(w1), w2_a(w2);
      const bool conj_a = helper.extendWIndices(w1_a, w2_a);
      const int i_a = b1 + nb * k1 + no * w1_a;
      const int j_a = b4 + nb * k2 + no * w2_a;

      int w1_b(helper.addWex(w2, w_ex));
      int w2_b(helper.addWex(w1, w_ex));
      const bool conj_b = helper.extendWIndices(w1_b, w2_b);
      const int i_b = b2 + nb * helper.addKex(k2, k_ex) + no * w1_b;
      const int j_b = b3 + nb * helper.addKex(k1, k_ex) + no * w2_b;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
      const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2) * (-factor);
    }
      {
        int w1_a(w1);
        int w2_a(helper.addWex(w1, w_ex));
        const bool conj_a = helper.extendWIndices(w1_a, w2_a);
        const int i_a = b1 + nb * k1 + no * w1_a;
        const int j_a = b3 + nb * helper.addKex(k1, k_ex) + no * w2_a;

        int w1_b(helper.addWex(w2, w_ex));
        int w2_b(w2);
        const bool conj_b = helper.extendWIndices(w1_b, w2_b);
        const int i_b = b2 + nb * helper.addKex(k2, k_ex) + no * w1_b;
        const int j_b = b4 + nb * k2 + no * w2_b;

        const CudaComplex<Real> Ga =
            cond_conj(G_up[i_a + ldgu * j_a] - G_down[i_a + ldgd * j_a], conj_a);
        const CudaComplex<Real> Gb =
            cond_conj(G_up[i_b + ldgu * j_b] - G_down[i_b + ldgd * j_b], conj_b);

        contribution += (Ga * Gb) * factor;
      }
      break;

    case PARTICLE_HOLE_CHARGE: {
      int w1_a(w1), w2_a(w2);
      const bool conj_a = helper.extendWIndices(w1_a, w2_a);
      const int i_a = b1 + nb * k1 + no * w1_a;
      const int j_a = b4 + nb * k2 + no * w2_a;

      int w1_b(helper.addWex(w2, w_ex));
      int w2_b(helper.addWex(w1, w_ex));
      const bool conj_b = helper.extendWIndices(w1_b, w2_b);
      const int i_b = b2 + nb * helper.addKex(k2, k_ex) + no * w1_b;
      const int j_b = b3 + nb * helper.addKex(k1, k_ex) + no * w2_b;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Gb_1 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
      const CudaComplex<Real> Gb_2 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2) * (-factor);
    }
      {
        int w1_a(w1);
        int w2_a(helper.addWex(w1, w_ex));
        const bool conj_a = helper.extendWIndices(w1_a, w2_a);
        const int i_a = b1 + nb * k1 + no * w1_a;
        const int j_a = b3 + nb * helper.addKex(k1, k_ex) + no * w2_a;

        int w1_b(helper.addWex(w2, w_ex));
        int w2_b(w2);
        const bool conj_b = helper.extendWIndices(w1_b, w2_b);
        const int i_b = b2 + nb * helper.addKex(k2, k_ex) + no * w1_b;
        const int j_b = b4 + nb * k2 + no * w2_b;

        const CudaComplex<Real> Ga =
            cond_conj(G_up[i_a + ldgu * j_a] + G_down[i_a + ldgd * j_a], conj_a);
        const CudaComplex<Real> Gb =
            cond_conj(G_up[i_b + ldgu * j_b] + G_down[i_b + ldgd * j_b], conj_b);

        contribution += (Ga * Gb) * factor;
      }
      break;

    case PARTICLE_PARTICLE_UP_DOWN: {
      int w1_a(w1), w2_a(w2);
      const bool conj_a = helper.extendWIndices(w1_a, w2_a);
      const int i_a = b1 + nb * k1 + no * w1_a;
      const int j_a = b3 + nb * k2 + no * w2_a;

      int w1_b(helper.wexMinus(w1, w_ex));
      int w2_b(helper.wexMinus(w2, w_ex));
      const bool conj_b = helper.extendWIndices(w1_b, w2_b);
      const int i_b = b2 + nb * helper.kexMinus(k1, k_ex) + no * w1_b;
      const int j_b = b4 + nb * helper.kexMinus(k2, k_ex) + no * w2_b;

      const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
      const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

      const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
      const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

      contribution = (Ga_1 * Gb_1 + Ga_2 * Gb_2) * factor;
    } break;
    default:  // abort
      asm("trap;");
  }
  CudaComplex<Real>* const result_ptr =
      G4 + helper.g4Index(b1, b2, b3, b4, k1, k2, k_ex, w1, w2, w_ex);
  dca::linalg::atomicAdd(result_ptr, contribution);
}

template <typename Real, FourPointType type>
void updateG4(std::complex<Real>* G4, const std::complex<Real>* G_up, const int ldgu,
              const std::complex<Real>* G_down, const int ldgd, const int nb, const int nk,
              const int nw_pos, const int nw_exchange, const int nk_exchange, const int sign,
              cudaStream_t stream) {
  const static int nw = 2 * nw_pos;
  const static int size_12 = nw * nk * nb * nb;
  const static int size_3 = nw_exchange * nk_exchange;
  const static auto blocks = getBlockSize3D(size_12, size_12, size_3);

  updateG4Kernel<Real, type><<<blocks[0], blocks[1], 0, stream>>>(
      castCudaComplex(G4), castCudaComplex(G_up), ldgu, castCudaComplex(G_down), ldgd, nb, nk, nw,
      nw_exchange, nk_exchange, sign, global::helper);

  // Check for errors.
  auto err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    linalg::util::printErrorMessage(err, __FUNCTION__, __FILE__, __LINE__);
    throw(std::runtime_error("CUDA failed to launch the G4 kernel."));
  }
}

void initializeG4Helpers(int nb, int nk, int nw_pos, const std::vector<int>& delta_k_indices,
                         const std::vector<int>& delta_w_indices, const int* add_k, int lda,
                         const int* sub_k, int lds) {
  if (!global::helper.isInitialized())
    global::helper.set(nb, nk, nw_pos, delta_k_indices, delta_w_indices, add_k, lda, sub_k, lds);
  assert(cudaPeekAtLastError() == cudaSuccess);
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

template void updateG4<float, PARTICLE_HOLE_TRANSVERSE>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);
template void updateG4<float, PARTICLE_HOLE_MAGNETIC>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);
template void updateG4<float, PARTICLE_HOLE_CHARGE>(std::complex<float>* G4,
                                                    const std::complex<float>* G_up, const int ldgu,
                                                    const std::complex<float>* G_down, const int ldgd,
                                                    const int nb, const int nk, const int nw_pos,
                                                    const int nw_exchange, const int nk_exchange,
                                                    const int sign, cudaStream_t stream);
template void updateG4<float, PARTICLE_PARTICLE_UP_DOWN>(
    std::complex<float>* G4, const std::complex<float>* G_up, const int ldgu,
    const std::complex<float>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);

template void updateG4<double, PARTICLE_HOLE_TRANSVERSE>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);
template void updateG4<double, PARTICLE_HOLE_MAGNETIC>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);
template void updateG4<double, PARTICLE_HOLE_CHARGE>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);
template void updateG4<double, PARTICLE_PARTICLE_UP_DOWN>(
    std::complex<double>* G4, const std::complex<double>* G_up, const int ldgu,
    const std::complex<double>* G_down, const int ldgd, const int nb, const int nk, const int nw_pos,
    const int nw_exchange, const int nk_exchange, const int sign, cudaStream_t stream);

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca
