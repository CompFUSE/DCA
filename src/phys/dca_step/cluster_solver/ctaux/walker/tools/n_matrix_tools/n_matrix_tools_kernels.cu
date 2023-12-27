// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements n_matrix_tools_kernels.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/n_matrix_tools/n_matrix_tools_kernels.hpp"

#include <cassert>

#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/gpu_type_mapping.hpp"
#include "dca/util/type_help.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace nkernels {
// dca::phys::solver::ctaux::nkernels::

using namespace dca::linalg;
using dca::util::RealAlias;
using dca::util::castGPUType;

const static int BLOCK_SIZE_x = 32;
const static int BLOCK_SIZE_y = 8;

template <class T>
__global__ void compute_G_cols_kernel(int N_i, int N_r, int N_c, const int* p_ptr,
                                      const T* exp_V_ptr, const T* N_ptr, int N_ld, const T* G_ptr,
                                      int G_ld, T* G_cols_ptr, int G_cols_ld) {
  int I = threadIdx.x + blockIdx.x * BLOCK_SIZE_x;  // blockDim.x;

  int l_MIN = BLOCK_SIZE_y * (blockIdx.y + 0);
  int l_MAX = BLOCK_SIZE_y * (blockIdx.y + 1);

  l_MIN = max(l_MIN, 0);
  l_MAX = min(l_MAX, N_i);

  if (I < N_r) {
    T the_one{};
    the_one += RealAlias<T>{1.0};
    for (int l = l_MIN; l < l_MAX; ++l) {
      if (p_ptr[l] >= N_c) {
        G_cols_ptr[I + l * G_cols_ld] = G_ptr[I + (p_ptr[l] - N_c) * G_ld];
      }
      else {
        const auto alpha = exp_V_ptr[l] / (exp_V_ptr[l] - the_one);

        G_cols_ptr[I + l * G_cols_ld] = alpha * N_ptr[I + p_ptr[l] * N_ld];
      }
    }

    // for(int l=0; l<N_i; ++l)
    for (int l = l_MIN; l < l_MAX; ++l)
      if (p_ptr[l] < N_c and I == p_ptr[l])
        G_cols_ptr[I + l * G_cols_ld] -= 1. / (exp_V_ptr[l] - the_one);
  }
}

template <class T>
void compute_G_cols(int N_i, int N_r, int N_c, const int* p_ptr, const T* exp_V_ptr, const T* N_ptr,
                    int N_ld, const T* G_ptr, int G_ld, T* G_cols_ptr, int G_cols_ld, int thread_id,
                    int stream_id) {
  if (N_r > 0 and N_i > 0) {
    int bl_x = dca::util::ceilDiv(N_r, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(N_i, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_G_cols_kernel<<<blocks, threads, 0, stream_handle>>>(
        N_i, N_r, N_c, p_ptr, castGPUType(exp_V_ptr), castGPUType(N_ptr), N_ld, castGPUType(G_ptr),
        G_ld, castGPUType(G_cols_ptr), G_cols_ld);

    checkErrorsCudaDebug();
  }
}

template void compute_G_cols(int, int, int, const int*, const float*, const float*, int,
                             const float*, int, float*, int, int, int);
template void compute_G_cols(int, int, int, const int*, const double*, const double*, int,
                             const double*, int, double*, int, int, int);
template void compute_G_cols(int, int, int, const int*, const std::complex<float>*,
                             const std::complex<float>*, int, const std::complex<float>*, int,
                             std::complex<float>*, int, int, int);
template void compute_G_cols(int, int, int, const int*, const std::complex<double>*,
                             const std::complex<double>*, int, const std::complex<double>*, int,
                             std::complex<double>*, int, int, int);

template <class T>
__global__ void compute_d_vector_kernel(int N_i, const int* d_ind, T* d_ptr, int* p_ptr,
                                        const T* N_ptr, int N_ld) {
  int I = threadIdx.x + blockIdx.x * blockDim.x;

  if (I < N_i) {
    int index = p_ptr[d_ind[I]];

    d_ptr[d_ind[I]] = 1. / N_ptr[index + index * N_ld];
  }
}

template <class T>
void compute_d_vector(int N_i, const int* d_ind, T* d_ptr, int* p_ptr, const T* N_ptr, int N_ld,
                      int thread_id, int stream_id) {
  if (N_i > 0) {
    checkErrorsCudaDebug();

    int Nr_t = 32;
    int Nr_b = dca::util::ceilDiv(N_i, Nr_t);

    dim3 threads(Nr_t);
    dim3 blocks(Nr_b);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    compute_d_vector_kernel<<<blocks, threads, 0, stream_handle>>>(N_i, d_ind, d_ptr, p_ptr, N_ptr,
                                                                   N_ld);

    checkErrorsCudaDebug();
  }
}

}  // namespace nkernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
