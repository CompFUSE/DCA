// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the GPU kernels used by the Dnfft1D class.

#include "dca/math/nfft/kernels_interface.hpp"

#include <array>
#include <cuda_runtime.h>

#include "dca/math/nfft/nfft_helper.cuh"
#include "dca/linalg/util/atomic_add_cuda.cu.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

std::array<int, 2> getBlockSize(const int ni, const int block_size) {
  const int n_threads = std::min(block_size, ni);
  const int n_blocks = util::ceilDiv(ni, n_threads);
  return std::array<int, 2>{n_blocks, n_threads};
}

std::array<dim3, 2> getBlockSize(const uint i, const uint j, const uint block_size = 32) {
  const uint n_threads_i = std::min(block_size, i);
  const uint n_threads_j = std::min(block_size, j);
  if (n_threads_i * n_threads_j > 32 * 32)
    throw(std::logic_error("Block size is too big"));

  const uint n_blocks_i = dca::util::ceilDiv(i, n_threads_i);
  const uint n_blocks_j = dca::util::ceilDiv(j, n_threads_j);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j), dim3(n_threads_i, n_threads_j)};
}

// TODO: consider constant or texture memory for the coefficients.
template <typename ScalarType, int oversampling, int window_sampling, bool accumulate_m_sqr = false>
__global__ void accumulateOnDeviceKernel(const ScalarType* M, const int ldm, const ScalarType sign,
                                         ScalarType* out, ScalarType* out_sqr, int ldo,
                                         const ConfigElem* config_left,
                                         const ConfigElem* config_right, const ScalarType* times,
                                         const ScalarType* cubic_coeff, const int size) {
  //  const int conv_size = 2 * oversampling + 1;
  const int thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_idx >= size * size)
    return;

  // Compute input indices.
  //  const int m_idx = thread_idx / conv_size;
  const int id_j = thread_idx / size;
  const int id_i = thread_idx - size * id_j;
  const ScalarType tau = nfft_helper.computeTau(times[id_i], times[id_j]);
  // Compute index of the convolution output index relative to a single input value.
  //  const int conv_idx = thread_idx - conv_size * m_idx;

  int t_idx, conv_coeff_idx;
  ScalarType delta_t;
  nfft_helper.computeInterpolationIndices<CUBIC, oversampling, window_sampling>(
      tau, t_idx, conv_coeff_idx, delta_t);

  const int linindex =
      nfft_helper.computeLinearIndex(config_left[id_i].band, config_right[id_j].band,
                                     config_left[id_i].site, config_right[id_j].site);
  ScalarType* const out_ptr = out + t_idx + ldo * linindex;

  const ScalarType f_val = sign * M[id_i + ldm * id_j];

  const ScalarType* conv_coeff = cubic_coeff + conv_coeff_idx;
  for (int l = 0; l < 2 * oversampling + 1; ++l) {
    const ScalarType conv_function_value =
        ((conv_coeff[3] * delta_t + conv_coeff[2]) * delta_t + conv_coeff[1]) * delta_t +
        conv_coeff[0];
    const ScalarType contribution = f_val * conv_function_value;
    linalg::atomicAdd(out_ptr, contribution);
    if (accumulate_m_sqr) {
      linalg::atomicAdd(out_sqr, f_val * conv_function_value * M[id_i + ldm * id_j]);
    }
    conv_coeff += 4;
  }
}

template <typename ScalarType, int oversampling, int window_sampling>
void accumulateOnDevice(const ScalarType* M, const int ldm, const ScalarType sign, ScalarType* out,
                        ScalarType* out_sqr, const int ldo, const ConfigElem* config_left,
                        const ConfigElem* config_right, const ScalarType* tau,
                        const ScalarType* cubic_coeff, const int size, cudaStream_t stream_) {
  const auto blocks = getBlockSize(size * size, 128);

  // TODO: check if there is a performance gain in using a block size that is a multiple of
  //       convolution_size.
  accumulateOnDeviceKernel<ScalarType, oversampling, window_sampling>
      <<<blocks[0], blocks[1], 0, stream_>>>(M, ldm, sign, out, out_sqr, ldo, config_left,
                                             config_right, tau, cubic_coeff, size);
}

template <typename ScalarType>
__global__ void sumKernel(const ScalarType* in, const int ldi, ScalarType* out, const int ldo,
                          const int n, const int m) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (i < n && j < m) {
    out[i + ldo * j] += in[i + ldi * j];
  }
}

template <typename ScalarType>
void sum(const ScalarType* in, const int ldi, ScalarType* out, const int ldo, const int n,
         const int m, cudaStream_t stream) {
  auto blocks = getBlockSize(n, m, 16);
  sumKernel<<<blocks[0], blocks[1], 0, stream>>>(in, ldi, out, ldo, n, m);
}

void initializeNfftHelper(int nb, int nc, const int* add_r, int lda, const int* sub_r, int lds,
                          double t0, double delta_t, double t0_window, double delta_t_window,
                          double beta) {
  NfftHelper::set(nb, nc, add_r, lda, sub_r, lds, t0, delta_t, t0_window, delta_t_window, beta);
}

// Explicit instantiation.
constexpr int oversampling = 8;
constexpr int window_sampling = 32;
template void accumulateOnDevice<double, oversampling, window_sampling>(
    const double* M, const int ldm, const double sign, double* out, double* out_sqr, const int ldo,
    const ConfigElem* config_left, const ConfigElem* config_right, const double* tau,
    const double* cubic_coeff, const int size, cudaStream_t stream_);
template void accumulateOnDevice<float, oversampling, window_sampling>(
    const float* M, const int ldm, const float sign, float* out, float* out_sqr, const int ldo,
    const ConfigElem* config_left, const ConfigElem* config_right, const float* tau,
    const float* cubic_coeff, const int size, cudaStream_t stream_);

template void sum<double>(const double* in, const int ldi, double* out, const int ldo, const int n,
                          const int m, cudaStream_t stream);
template void sum<float>(const float* in, const int ldi, float* out, const int ldo, const int n,
                         const int m, cudaStream_t stream);

}  // namespace details
}  // namespace nfft
}  // namespace math
}  // namespace dca
