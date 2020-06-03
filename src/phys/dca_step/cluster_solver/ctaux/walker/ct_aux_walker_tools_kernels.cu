// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements ct_aux_walker_tools_kernels.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/ct_aux_walker_tools_kernels.hpp"

#include "cuda_runtime.h"

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace walkerkernels {
// dca::phys::solver::ctaux::walkerkernels::

template <typename Real>
__global__ void compute_Gamma_kernel(Real* Gamma, int Gamma_n, int Gamma_ld, Real* N, int N_r,
                                     int N_c, int N_ld, Real* G, int G_r, int G_c, int G_ld,
                                     int* random_vertex_vector, Real* exp_V, Real* exp_delta_V) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int vertex_index = N_c - G_c;

  if (i < Gamma_n and j < Gamma_n) {
    int configuration_e_spin_index_i = random_vertex_vector[i];
    int configuration_e_spin_index_j = random_vertex_vector[j];

    if (configuration_e_spin_index_j < vertex_index) {
      Real delta = 0;

      if (configuration_e_spin_index_i == configuration_e_spin_index_j)
        delta = 1.;

      Real N_ij = N[configuration_e_spin_index_i + configuration_e_spin_index_j * N_ld];

      Gamma[i + j * Gamma_ld] = (N_ij * exp_V[j] - delta) / (exp_V[j] - 1.);
    }
    else
      Gamma[i + j * Gamma_ld] =
          G[configuration_e_spin_index_i + (configuration_e_spin_index_j - vertex_index) * G_ld];
  }

  if (i < Gamma_n and j < Gamma_n and i == j) {
    Real gamma_k = exp_delta_V[j];
    Gamma[i + j * Gamma_ld] -= (gamma_k) / (gamma_k - 1.);
  }
}

template <typename Real>
void compute_Gamma(Real* Gamma, int Gamma_n, int Gamma_ld, Real* N, int N_r, int N_c, int N_ld,
                   Real* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, Real* exp_V,
                   Real* exp_delta_V, int thread_id, int stream_id) {
  const int number_of_threads = 16;

  if (Gamma_n > 0) {
    checkErrorsCudaDebug();

    dim3 threads(number_of_threads, number_of_threads);

    dim3 blocks(dca::util::ceilDiv(Gamma_n, number_of_threads),
                dca::util::ceilDiv(Gamma_n, number_of_threads));

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_Gamma_kernel<<<blocks, threads, 0, stream_handle>>>(
        Gamma, Gamma_n, Gamma_ld, N, N_r, N_c, N_ld, G, G_r, G_c, G_ld, random_vertex_vector, exp_V,
        exp_delta_V);

    checkErrorsCudaDebug();
  }
}
template void compute_Gamma(float* Gamma, int Gamma_n, int Gamma_ld, float* N, int N_r, int N_c,
                            int N_ld, float* G, int G_r, int G_c, int G_ld, int* random_vertex_vector,
                            float* exp_V, float* exp_delta_V, int thread_id, int stream_id);
template void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c,
                            int N_ld, double* G, int G_r, int G_c, int G_ld,
                            int* random_vertex_vector, double* exp_V, double* exp_delta_V,
                            int thread_id, int stream_id);

}  // namespace walkerkernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
