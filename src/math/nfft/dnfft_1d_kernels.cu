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

#include "dca/math/nfft/nfft_helper.cu.hpp"
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

template <typename ScalarType>
struct HelperSelector {
  static NfftHelperManager<ScalarType> value;
};
template <typename ScalarType>
NfftHelperManager<ScalarType> HelperSelector<ScalarType>::value;

template <typename ScalarType>
__global__ void accumulateOnDeviceKernel(const double* M, const int ldm, const int sign,
                                         ScalarType* out, ScalarType* out_sqr, int ldo,
                                         const ConfigElem* config_left,
                                         const ConfigElem* config_right, const ScalarType* times,
                                         const ScalarType* cubic_coeff, const int size,
                                         const NfftHelper<ScalarType> helper) {
  const int conv_size = 2 * helper.get_oversampling() + 1;
  const int thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_idx >= size * size * conv_size)
    return;

  // Compute input indices.
  const int m_idx = thread_idx / conv_size;
  const int id_j = m_idx / size;
  const int id_i = m_idx - size * id_j;
  const ScalarType f_val = M[id_i + ldm * id_j];
  const int linindex = helper.computeLinearIndex(config_left[id_i].band, config_right[id_j].band,
                                                 config_left[id_i].site, config_right[id_j].site);
  const ScalarType tau = helper.computeTau(times[id_i], times[id_j]);
  // Compute index of the convolution output index relative to a single input value.
  const int conv_idx = thread_idx - conv_size * m_idx;

  int t_idx, conv_coeff_idx;
  ScalarType delta_t;
  helper.computeInterpolationIndices<CUBIC>(tau, t_idx, conv_coeff_idx, delta_t);

  const ScalarType* conv_coeff = cubic_coeff + conv_coeff_idx + 4 * conv_idx;
  const ScalarType conv_function_value = conv_coeff[0] + conv_coeff[1] * delta_t +
                                         conv_coeff[2] * delta_t * delta_t +
                                         conv_coeff[3] * delta_t * delta_t * delta_t;

  ScalarType* const out_ptr = out + t_idx + conv_idx + ldo * linindex;
  linalg::atomicAdd(out_ptr, sign * f_val * conv_function_value);

  if (out_sqr != nullptr) {
    ScalarType* const out_sqr_ptr = out_sqr + (out_ptr - out);
    linalg::atomicAdd(out_sqr_ptr, sign * f_val * f_val * conv_function_value);
  }
}

template <typename ScalarType>
void accumulateOnDevice(const double* M, const int ldm, const int sign, ScalarType* out,
                        ScalarType* out_sqr, const int ldo, const ConfigElem* config_left,
                        const ConfigElem* config_right, const ScalarType* tau,
                        const ScalarType* cubic_coeff, const int size, cudaStream_t stream_) {
  const auto& helper = HelperSelector<ScalarType>::value;
  const static int convolution_size = 2 * helper.get_oversampling() + 1;
  const auto blocks = getBlockSize(size * size * convolution_size, 85);

  // TODO: check if there is a performance gain in using a block size that is a multiple of
  //       convolution_size.
  accumulateOnDeviceKernel<ScalarType><<<blocks[0], blocks[1], 0, stream_>>>(
      M, ldm, sign, out, out_sqr, ldo, config_left, config_right, tau, cubic_coeff, size, helper);
}

template <typename ScalarType>
__global__ void sumKernel(const ScalarType* in, const int ldi, ScalarType* out, const int ldo,
                          const int n, const int m) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  out[i + ldo * j] += in[i + ldi * j];
}

template <typename ScalarType>
void sum(const ScalarType* in, const int ldi, ScalarType* out, const int ldo, const int n,
         const int m, cudaStream_t stream) {
  auto blocks = getBlockSize(n, m, 16);
  sumKernel<<<blocks[0], blocks[1], 0, stream>>>(in, ldi, out, ldo, n, m);
}

template <typename ScalarType>
void initializeNfftHelper(int nb, int nr, const int* sub_r, int lds, int oversampling,
                          int window_sampling, ScalarType t0, ScalarType delta_t,
                          ScalarType t0_window, ScalarType delta_t_window, ScalarType beta) {
  auto& helper = HelperSelector<ScalarType>::value;
  if (helper.isInitialized())
    return;

  helper.set(nb, nr, sub_r, lds, oversampling, window_sampling, t0, delta_t, t0_window,
             delta_t_window, beta);
}

// Explicit instantiation.
template void accumulateOnDevice<double>(const double* M, const int ldm, const int sign,
                                         double* out, double* out_sqr, const int ldo,
                                         const ConfigElem* config_left,
                                         const ConfigElem* config_right, const double* tau,
                                         const double* cubic_coeff, const int size,
                                         cudaStream_t stream_);
template void accumulateOnDevice<float>(const double* M, const int ldm, const int sign, float* out,
                                        float* out_sqr, const int ldo, const ConfigElem* config_left,
                                        const ConfigElem* config_right, const float* tau,
                                        const float* cubic_coeff, const int size,
                                        cudaStream_t stream_);
template void sum<double>(const double* in, const int ldi, double* out, const int ldo, const int n,
                          const int m, cudaStream_t stream);
template void sum<float>(const float* in, const int ldi, float* out, const int ldo, const int n,
                         const int m, cudaStream_t stream);
template void initializeNfftHelper<double>(int nb, int nr, const int* sub_r, int lds,
                                           int oversampling, int window_sampling, double t0,
                                           double delta_t, double t0_window, double delta_t_window,
                                           double beta);
template void initializeNfftHelper<float>(int nb, int nr, const int* sub_r, int lds, int oversampling,
                                          int window_sampling, float t0, float delta_t,
                                          float t0_window, float delta_t_window, float beta);

}  // details
}  // nfft
}  // math
}  // dca
