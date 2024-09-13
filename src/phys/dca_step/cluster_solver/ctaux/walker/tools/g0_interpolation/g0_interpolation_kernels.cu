// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements g0_interpolation_kernels.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_kernels.hpp"

#include <cassert>
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_types.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace g0kernels {
// dca::phys::solver::ctaux::g0kernels::

__global__ void interpolation_kernel(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                     double* t, double* G0, std::pair<int, int> G0_cs,
                                     std::pair<int, int> G0_gs, double* r0_min_r1,
                                     std::pair<int, int> r0_min_r1_cs,
                                     std::pair<int, int> r0_min_r1_gs, double* G0_r_t,
                                     std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                     double* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                     std::pair<int, int> grad_G0_r_t_gs) {
  int I = threadIdx.x + blockDim.x * blockIdx.x;
  int J = threadIdx.y + blockDim.y * blockIdx.y;

  if (I > -1 && I < Nv && J > -1 && J < Nv) {
    int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
    double tau = t[I] - t[J];

    double scaled_tau = (tau + beta) * double(Nt) / (2. * beta);

    int t_ind = scaled_tau;
    double delta_tau = scaled_tau - t_ind;

    int linind = b[I] + Nb * (b[J] + Nb * delta_r);

    double f_0 = G0_r_t[t_ind + linind * G0_r_t_gs.first];
    double grad = grad_G0_r_t[t_ind + linind * grad_G0_r_t_gs.first];

    G0[I + G0_gs.first * J] = -(f_0 + grad * delta_tau);

    if (I == J) {
      int t0_ind = Nt / 2;
      int r0_ind = r0_min_r1[0];

      G0[I + G0_gs.first * I] =
          -G0_r_t[t0_ind + G0_r_t_gs.first * (b[I] + Nb * (b[I] + Nb * r0_ind))];
    }
  }
}

void interpolate_G0_matrix_on_GPU(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                  double* t, double* G0, std::pair<int, int> G0_cs,
                                  std::pair<int, int> G0_gs, double* r0_min_r1,
                                  std::pair<int, int> r0_min_r1_cs,
                                  std::pair<int, int> r0_min_r1_gs, double* G0_r_t,
                                  std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                  double* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                  std::pair<int, int> grad_G0_r_t_gs) {
  checkErrorsCudaDebug();

  int Nr_t = 16;
  int Nr_b = dca::util::ceilDiv(Nv, Nr_t);

  dim3 threads(Nr_t, Nr_t);
  dim3 blocks(Nr_b, Nr_b);

  interpolation_kernel<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nv, b, r, t, G0, G0_cs, G0_gs,
                                            r0_min_r1, r0_min_r1_cs, r0_min_r1_gs, G0_r_t, G0_r_t_cs,
                                            G0_r_t_gs, grad_G0_r_t, grad_G0_r_t_cs, grad_G0_r_t_gs);

  checkErrorsCudaDebug();
}


/***********************************
 ***
 ***        AKIMA-interpolation
 ***
 ************************************/

using namespace dca::linalg;
using dca::util::castGPUType;

const static int BLOCK_SIZE_x = 32;
const static int BLOCK_SIZE_y = 16;

template <typename Scalar, typename Real>
__global__ void akima_interpolation_fat_column(
    int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, const int* b, const int* r, const Real* t,
    Scalar* G0, std::pair<int, int> G0_cs, std::pair<int, int> G0_gs, const int* r0_min_r1,
    std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs, const Scalar* alpha,
    std::pair<int, int> alpha_cs, std::pair<int, int> alpha_gs) {
  assert(blockDim.x == BLOCK_SIZE_x);

  int I = threadIdx.x + BLOCK_SIZE_x * blockIdx.x;

  int J_min = Nc + BLOCK_SIZE_y * (blockIdx.y + 0);
  int J_max = Nc + BLOCK_SIZE_y * (blockIdx.y + 1);

  J_min = max(Nc, J_min);
  J_max = min(Nv, J_max);

  if (I > -1 && I < Nv) {
    for (int J = J_min; J < J_max; ++J) {
      int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
      Real tau = t[I] - t[J];

      Real scaled_tau = (tau + beta) * Real(Nt) / (2. * beta);

      int t_ind = scaled_tau;
      Real delta_tau = scaled_tau - t_ind;

      assert(delta_tau > -1.e-16 and delta_tau < 1 + -1.e-16);

      int col_ind = b[I] + Nb * (b[J] + Nb * delta_r);
      int row_ind = 4 * t_ind;
      int alpha_LD = alpha_gs.first;

      assert(row_ind > -1 and row_ind < alpha_cs.first);
      assert(col_ind > -1 and col_ind < alpha_cs.second);

      const Scalar* a_ptr = &alpha[row_ind + col_ind * alpha_LD];

      assert(I > -1 and I < G0_cs.first);
      assert(J > -1 and J < G0_cs.second);

      G0[I + G0_gs.first * J] =
          -(a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));
    }
  }

  if (I >= J_min && I < J_max)  // I==J
  {
    int t0_ind = Nt / 2;
    int r0_ind = r0_min_r1[0];

    int col_ind = b[I] + Nb * (b[I] + Nb * r0_ind);
    int row_ind = 4 * t0_ind;
    int alpha_LD = alpha_gs.first;

    assert(row_ind > -1 and row_ind < alpha_cs.first);
    assert(col_ind > -1 and col_ind < alpha_cs.second);

    assert(I > -1 and I < G0_cs.first);

    G0[I + G0_gs.first * I] = -alpha[row_ind + col_ind * alpha_LD];
  }
}

template <typename Scalar, typename Real>
__global__ void akima_interpolation_fat_row(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv,
                                            const int* b, const int* r, const Real* t, Scalar* G0,
                                            std::pair<int, int> G0_cs, std::pair<int, int> G0_gs,
                                            const int* r0_min_r1, std::pair<int, int> r0_min_r1_cs,
                                            std::pair<int, int> r0_min_r1_gs, const Scalar* alpha,
                                            std::pair<int, int> alpha_cs,
                                            std::pair<int, int> alpha_gs) {
  static_assert(std::is_same<Real, float>::value || std::is_same<Real, double>::value, "");

  int I = Nc + threadIdx.x + BLOCK_SIZE_x * blockIdx.x;

  int J_min = BLOCK_SIZE_y * (blockIdx.y + 0);
  int J_max = BLOCK_SIZE_y * (blockIdx.y + 1);

  J_min = max(J_min, 0);
  J_max = min(J_max, Nc);

  if (I >= Nc && I < Nv) {
    for (int J = J_min; J < J_max; ++J) {
      int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
      Real tau = t[I] - t[J];

      Real scaled_tau = (tau + beta) * Real(Nt) / (2. * beta);

      const int t_ind = scaled_tau;
      const Real delta_tau = scaled_tau - t_ind;
      assert(delta_tau > -1.e-16 && delta_tau < 1 + -1.e-16);

      int col_ind = b[I] + Nb * (b[J] + Nb * delta_r);
      int row_ind = 4 * t_ind;
      int alpha_LD = alpha_gs.first;

      assert(row_ind > -1 && row_ind < alpha_cs.first);
      assert(col_ind > -1 && col_ind < alpha_cs.second);

      const Scalar* const a_ptr = &alpha[row_ind + col_ind * alpha_LD];

      assert(I > -1 && I < G0_cs.first);
      assert(J > -1 && J < G0_cs.second);

      G0[I + G0_gs.first * J] =
          -(a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));
    }
  }
}

template <typename Scalar, typename Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, const int* b,
                                const int* r, const Real* t, Scalar* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, const int* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                const Scalar* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs) {
  akima_interpolation_on_GPU(Nb, Nr, Nt, beta, Nc, Nv, b, r, t, G0, G0_cs, G0_gs, r0_min_r1,
                             r0_min_r1_cs, r0_min_r1_gs, alpha, alpha_cs, alpha_gs, 0, 0);

  dca::linalg::util::getStream(0, 0).sync();
}

// Note: two template parameters are used as the nvcc linker does not seem to behave correctly when
// instantiating a template with dependent types.
template <typename Scalar, typename Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, const int* b,
                                const int* r, const Real* t, Scalar* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, const int* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                const Scalar* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs, int thread_id, int stream_id) {
  static_assert(std::is_same<dca::util::RealAlias<Scalar>, Real>::value,
                "Real is not the real type of Scalar");
  if (Nv - Nc > 0 && Nv > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    akima_interpolation_fat_column<<<blocks, threads, 0, stream_handle>>>(
        Nb, Nr, Nt, beta, Nc, Nv, b, r, t, castGPUType(G0), G0_cs, G0_gs, r0_min_r1, r0_min_r1_cs,
        r0_min_r1_gs, castGPUType(alpha), alpha_cs, alpha_gs);

    checkErrorsCudaDebug();
  }

  if (Nv - Nc > 0 and Nc > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    akima_interpolation_fat_row<<<blocks, threads, 0, stream_handle>>>(
        Nb, Nr, Nt, beta, Nc, Nv, b, r, t, castGPUType(G0), G0_cs, G0_gs, r0_min_r1, r0_min_r1_cs,
        r0_min_r1_gs, castGPUType(alpha), alpha_cs, alpha_gs);

    checkErrorsCudaDebug();
  }
}

// Instantiation.
template void akima_interpolation_on_GPU(int, int, int, float, int, int, const int*, const int*,
                                         const float*, float*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const float*, std::pair<int, int>,
                                         std::pair<int, int>, int, int);
template void akima_interpolation_on_GPU(int, int, int, double, int, int, const int*, const int*,
                                         const double*, double*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const double*, std::pair<int, int>,
                                         std::pair<int, int>, int, int);
template void akima_interpolation_on_GPU(int, int, int, float, int, int, const int*, const int*,
                                         const float*, std::complex<float>*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const std::complex<float>*,
                                         std::pair<int, int>, std::pair<int, int>, int, int);
template void akima_interpolation_on_GPU(int, int, int, double, int, int, const int*, const int*,
                                         const double*, std::complex<double>*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const std::complex<double>*,
                                         std::pair<int, int>, std::pair<int, int>, int, int);

template void akima_interpolation_on_GPU(int, int, int, float, int, int, const int*, const int*,
                                         const float*, float*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const float*, std::pair<int, int>,
                                         std::pair<int, int>);
template void akima_interpolation_on_GPU(int, int, int, double, int, int, const int*, const int*,
                                         const double*, double*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const double*, std::pair<int, int>,
                                         std::pair<int, int>);
template void akima_interpolation_on_GPU(int, int, int, float, int, int, const int*, const int*,
                                         const float*, std::complex<float>*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const std::complex<float>*,
                                         std::pair<int, int>, std::pair<int, int>);
template void akima_interpolation_on_GPU(int, int, int, double, int, int, const int*, const int*,
                                         const double*, std::complex<double>*, std::pair<int, int>,
                                         std::pair<int, int>, const int*, std::pair<int, int>,
                                         std::pair<int, int>, const std::complex<double>*,
                                         std::pair<int, int>, std::pair<int, int>);

}  // namespace g0kernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
