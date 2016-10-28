// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for N-matrix tools.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_CU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_CU_HPP

#include <cassert>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace nkernels {
// dca::phys::solver::ctaux::nkernels::

const static int BLOCK_SIZE_x = 32;
const static int BLOCK_SIZE_y = 8;

__global__ void compute_G_cols_kernel(int N_i, int N_r, int N_c, int* p_ptr, double* exp_V_ptr,
                                      double* N_ptr, int N_ld, double* G_ptr, int G_ld,
                                      double* G_cols_ptr, int G_cols_ld) {
  int I = threadIdx.x + blockIdx.x * BLOCK_SIZE_x;  // blockDim.x;

  int l_MIN = BLOCK_SIZE_y * (blockIdx.y + 0);
  int l_MAX = BLOCK_SIZE_y * (blockIdx.y + 1);

  l_MIN = max(l_MIN, 0);
  l_MAX = min(l_MAX, N_i);

  if (I < N_r) {
    // for(int l=0; l<N_i; ++l)
    for (int l = l_MIN; l < l_MAX; ++l) {
      if (p_ptr[l] >= N_c) {
        G_cols_ptr[I + l * G_cols_ld] = G_ptr[I + (p_ptr[l] - N_c) * G_ld];
      }
      else {
        double alpha = exp_V_ptr[l] / (exp_V_ptr[l] - 1.);

        G_cols_ptr[I + l * G_cols_ld] = alpha * N_ptr[I + p_ptr[l] * N_ld];
      }
    }

    // for(int l=0; l<N_i; ++l)
    for (int l = l_MIN; l < l_MAX; ++l)
      if (p_ptr[l] < N_c and I == p_ptr[l])
        G_cols_ptr[I + l * G_cols_ld] -= 1. / (exp_V_ptr[l] - 1.);
  }
}

void compute_G_cols(int N_i, int N_r, int N_c, int* p_ptr, double* exp_V_ptr, double* N_ptr,
                    int N_ld, double* G_ptr, int G_ld, double* G_cols_ptr, int G_cols_ld,
                    int thread_id, int stream_id) {
  if (N_r > 0 and N_i > 0) {
    int bl_x = dca::util::ceilDiv(N_r, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(N_i, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_G_cols_kernel<<<blocks, threads, 0, stream_handle>>>(
        N_i, N_r, N_c, p_ptr, exp_V_ptr, N_ptr, N_ld, G_ptr, G_ld, G_cols_ptr, G_cols_ld);

#ifdef DEBUG_CUDA
    cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
  }
}

__global__ void compute_d_vector_kernel(int N_i, int* d_ind, double* d_ptr, int* p_ptr,
                                        double* N_ptr, int N_ld) {
  int I = threadIdx.x + blockIdx.x * blockDim.x;

  if (I < N_i) {
    int index = p_ptr[d_ind[I]];

    d_ptr[d_ind[I]] = 1. / N_ptr[index + index * N_ld];
  }
}

void compute_d_vector(int N_i, int* d_ind, double* d_ptr, int* p_ptr, double* N_ptr, int N_ld,
                      int thread_id, int stream_id) {
  if (N_i > 0) {
#ifdef DEBUG_CUDA
    cuda_check_for_errors_bgn(__FUNCTION__, __FILE__, __LINE__);
#endif

    int Nr_t = 32;
    int Nr_b = dca::util::ceilDiv(N_i, Nr_t);

    dim3 threads(Nr_t);
    dim3 blocks(Nr_b);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_d_vector_kernel<<<blocks, threads, 0, stream_handle>>>(N_i, d_ind, d_ptr, p_ptr, N_ptr,
                                                                   N_ld);

#ifdef DEBUG_CUDA
    cuda_check_for_errors_end(__FUNCTION__, __FILE__, __LINE__);
#endif
  }
}

}  // nkernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_CU_HPP
