// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_CU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_CU_H

namespace DCA {
namespace QMCI {
namespace CT_AUX_WALKER_GPU_KERNELS {
// DCA::QMCI::CT_AUX_WALKER_GPU_KERNELS::

__global__ void compute_Gamma_kernel(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r,
                                     int N_c, int N_ld, double* G, int G_r, int G_c, int G_ld,
                                     int* random_vertex_vector, double* exp_V, double* exp_delta_V) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  int vertex_index = N_c - G_c;

  if (i < Gamma_n and j < Gamma_n) {
    int configuration_e_spin_index_i = random_vertex_vector[i];
    int configuration_e_spin_index_j = random_vertex_vector[j];

    if (configuration_e_spin_index_j < vertex_index) {
      double delta = 0;

      if (configuration_e_spin_index_i == configuration_e_spin_index_j)
        delta = 1.;

      double N_ij = N[configuration_e_spin_index_i + configuration_e_spin_index_j * N_ld];

      Gamma[i + j * Gamma_ld] = (N_ij * exp_V[j] - delta) / (exp_V[j] - 1.);
    }
    else
      Gamma[i + j * Gamma_ld] =
          G[configuration_e_spin_index_i + (configuration_e_spin_index_j - vertex_index) * G_ld];
  }

  if (i < Gamma_n and j < Gamma_n and i == j) {
    double gamma_k = exp_delta_V[j];
    Gamma[i + j * Gamma_ld] -= (gamma_k) / (gamma_k - 1.);
  }
}

void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c, int N_ld,
                   double* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, double* exp_V,
                   double* exp_delta_V, int thread_id, int stream_id) {
  const int number_of_threads = 16;

  if (Gamma_n > 0) {
#ifdef DEBUG_CUDA
    cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

    dim3 threads(number_of_threads, number_of_threads);

    dim3 blocks(dca::util::ceilDiv(Gamma_n, number_of_threads),
                dca::util::ceilDiv(Gamma_n, number_of_threads));

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_Gamma_kernel<<<blocks, threads, 0, stream_handle>>>(
        Gamma, Gamma_n, Gamma_ld, N, N_r, N_c, N_ld, G, G_r, G_c, G_ld, random_vertex_vector, exp_V,
        exp_delta_V);

#ifdef DEBUG_CUDA
    cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
  }
}

}  // CT_AUX_WALKER_GPU_KERNELS
}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_CU_H
