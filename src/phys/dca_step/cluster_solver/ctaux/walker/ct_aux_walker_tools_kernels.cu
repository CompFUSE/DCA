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

#include <complex>

#include "cuda_runtime.h"

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace walkerkernels {
// dca::phys::solver::ctaux::walkerkernels::

template <class T>
__global__ void compute_Gamma_kernel(T* Gamma, int Gamma_n, int Gamma_ld, const T* N, int N_r,
                                     int N_c, int N_ld, const T* G, int G_r, int G_c, int G_ld,
                                     const int* random_vertex_vector, const T* exp_V,
                                     const T* exp_delta_V) {
  using namespace dca::linalg;
  using dca::util::Real;

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  const int vertex_index = N_c - G_c;

  if (i < Gamma_n and j < Gamma_n) {
    const int configuration_e_spin_index_i = random_vertex_vector[i];
    const int configuration_e_spin_index_j = random_vertex_vector[j];

    if (configuration_e_spin_index_j < vertex_index) {
      Real<T> delta = 0;

      if (configuration_e_spin_index_i == configuration_e_spin_index_j)
        delta = 1.;

      const auto N_ij = N[configuration_e_spin_index_i + configuration_e_spin_index_j * N_ld];

      Gamma[i + j * Gamma_ld] = (N_ij * exp_V[j] - delta) / (exp_V[j] - 1.);
    }
    else
      Gamma[i + j * Gamma_ld] =
          G[configuration_e_spin_index_i + (configuration_e_spin_index_j - vertex_index) * G_ld];
  }

  if (i < Gamma_n and j < Gamma_n and i == j) {
    const auto gamma_k = exp_delta_V[j];
    Gamma[i + j * Gamma_ld] -= (gamma_k) / (gamma_k - 1.);
  }
}

template <class T>
void compute_Gamma(T* Gamma, int Gamma_n, int Gamma_ld, const T* N, int N_r, int N_c, int N_ld,
                   const T* G, int G_r, int G_c, int G_ld, const int* random_vertex_vector,
                   const T* exp_V, const T* exp_delta_V, int thread_id, int stream_id) {
  const int number_of_threads = 16;

  if (Gamma_n > 0) {
    checkErrorsCudaDebug();

    dim3 threads(number_of_threads, number_of_threads);

    dim3 blocks(dca::util::ceilDiv(Gamma_n, number_of_threads),
                dca::util::ceilDiv(Gamma_n, number_of_threads));

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    using dca::linalg::util::castCudaComplex;
    compute_Gamma_kernel<<<blocks, threads, 0, stream_handle>>>(
        castCudaComplex(Gamma), Gamma_n, Gamma_ld, castCudaComplex(N), N_r, N_c, N_ld,
        castCudaComplex(G), G_r, G_c, G_ld, random_vertex_vector, castCudaComplex(exp_V),
        castCudaComplex(exp_delta_V));

    checkErrorsCudaDebug();
  }
}

template void compute_Gamma(float*, int, int, const float*, int, int, int, const float*, int, int,
                            int, const int*, const float*, const float*, int, int);
template void compute_Gamma(double*, int, int, const double*, int, int, int, const double*, int,
                            int, int, const int*, const double*, const double*, int, int);
template void compute_Gamma(std::complex<float>*, int, int, const std::complex<float>*, int, int,
                            int, const std::complex<float>*, int, int, int, const int*,
                            const std::complex<float>*, const std::complex<float>*, int, int);
template void compute_Gamma(std::complex<double>*, int, int, const std::complex<double>*, int, int,
                            int, const std::complex<double>*, int, int, int, const int*,
                            const std::complex<double>*, const std::complex<double>*, int, int);

}  // namespace walkerkernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
