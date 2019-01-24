// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Helper class for computing interpolation indices during the single particle accumulation.

#ifndef DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER
#define DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER

#include <array>
#include <cassert>
#include <cuda.h>
#include <stdexcept>

#include "dca/math/nfft/nfft_mode_names.hpp"

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

template <typename ScalarType>
class NfftHelper {
public:
  __host__ void set(int nb, int nr, const int* sub_r, int lds, int oversampling,
                    int window_sampling, ScalarType t0, ScalarType delta_t, ScalarType t0_window,
                    ScalarType delta_t_window, ScalarType beta);

  NfftHelper(const NfftHelper& other) = default;

  __host__ bool isInitialized() const {
    return parameters_int_ != nullptr;
  }

  __device__ int inline computeLinearIndex(int b1, int b2, int r1, int r2) const;

  __device__ ScalarType inline computeTau(ScalarType t1, ScalarType t2) const;

  __device__ inline int deltaR(int r1, int r2) const;

  template <NfftModeNames mode>
  __device__ inline void computeInterpolationIndices(ScalarType t, int& t_idx, int& conv_coeff_idx,
                                                     ScalarType& delta_t) const;

  __device__ inline int get_oversampling() const {
    return parameters_int_[3];
  }

protected:
  // This object can be constructed only through its derived class.
  NfftHelper() = default;

  int* sub_matrix_ = nullptr;

  // Stores in this order { nb, nr, lds, oversampling, window_sampling}
  int* parameters_int_ = nullptr;
  // Stores in this order { t0, delta_t, t0_window, delta_t_window, one_div_delta,
  // one_div_delta_window, tau_scale}
  ScalarType* parameters_real_ = nullptr;
};

template <typename ScalarType>
class NfftHelperManager : public NfftHelper<ScalarType> {
public:
  NfftHelperManager() = default;

  __host__ ~NfftHelperManager();

  __host__ __device__ NfftHelperManager(const NfftHelperManager& other) = delete;

  __host__ void set(int nb, int nr, const int* sub_r, int lds, int oversampling,
                    int window_sampling, ScalarType t0, ScalarType delta_t, ScalarType t0_window,
                    ScalarType delta_t_window, ScalarType beta);

  int get_oversampling() const {
    assert(NfftHelper<ScalarType>::isInitialized());
    return oversampling_;
  }

private:
  int oversampling_ = -1;
};

template <typename ScalarType>
__host__ void NfftHelper<ScalarType>::set(const int nb, const int nr, const int* sub_r, int lds,
                                          int oversampling, int window_sampling, ScalarType t0,
                                          ScalarType delta_t, const ScalarType t0_window,
                                          const ScalarType delta_t_window, const ScalarType beta) {
  if (isInitialized())
    throw(std::logic_error("already initialized."));

  cudaMalloc(&sub_matrix_, sizeof(int) * lds * nr);
  cudaMemcpy(sub_matrix_, sub_r, sizeof(int) * lds * nr, cudaMemcpyHostToDevice);

  const std::array<int, 5> parameters_int_host{nb, nr, lds, oversampling, window_sampling};
  cudaMalloc(&parameters_int_, sizeof(int) * parameters_int_host.size());
  cudaMemcpy(parameters_int_, parameters_int_host.data(), sizeof(int) * parameters_int_host.size(),
             cudaMemcpyHostToDevice);

  ScalarType one(1);
  const ScalarType tau_scale = ScalarType(0.5) / beta;
  const std::array<ScalarType, 7> parameters_real_host{
      t0, delta_t, t0_window, delta_t_window, one / delta_t, one / delta_t_window, tau_scale};
  cudaMalloc(&parameters_real_, sizeof(ScalarType) * parameters_real_host.size());
  cudaMemcpy(parameters_real_, parameters_real_host.data(),
             sizeof(ScalarType) * parameters_real_host.size(), cudaMemcpyHostToDevice);
}

template <typename ScalarType>
__device__ int NfftHelper<ScalarType>::deltaR(const int r1, const int r2) const {
  const int ld = parameters_int_[2];
  return sub_matrix_[r1 + ld * r2];
}

template <typename ScalarType>
__host__ void NfftHelperManager<ScalarType>::set(const int nb, const int nr, const int* sub_r,
                                                 int lds, int oversampling, int window_sampling,
                                                 ScalarType t0, ScalarType delta_t,
                                                 const ScalarType t0_window,
                                                 const ScalarType delta_t_window,
                                                 const ScalarType beta) {
  NfftHelper<ScalarType>::set(nb, nr, sub_r, lds, oversampling, window_sampling, t0, delta_t,
                              t0_window, delta_t_window, beta);
  oversampling_ = oversampling;
}

template <typename ScalarType>
__device__ int inline NfftHelper<ScalarType>::computeLinearIndex(const int b1, const int b2,
                                                                 const int r1, const int r2) const {
  const int delta_r = deltaR(r2, r1);
  const int nb = parameters_int_[0];
  return b1 + nb * (b2 + nb * delta_r);
}

template <typename ScalarType>
__device__ ScalarType NfftHelper<ScalarType>::computeTau(ScalarType t1, ScalarType t2) const {
  const ScalarType tau_scale = parameters_real_[6];
  return (t1 - t2) * tau_scale;
}

template <typename ScalarType>
template <NfftModeNames mode>
__device__ void NfftHelper<ScalarType>::computeInterpolationIndices(ScalarType t, int& t_idx,
                                                                    int& conv_coeff_idx,
                                                                    ScalarType& delta_t) const {
  static_assert(mode == CUBIC || mode == LINEAR,
                "This method is defined only for linear or cubic interpolations.");

  const ScalarType t0 = parameters_real_[0];
  const ScalarType delta_t_padded = parameters_real_[1];
  const ScalarType t0_window = parameters_real_[2];
  const ScalarType delta_t_window = parameters_real_[3];
  auto get_tau = [=](int index) { return t0 + index * delta_t_padded; };
  auto get_window_tau = [&](int index) { return t0_window + index * delta_t_window; };
  const ScalarType one_div_delta = parameters_real_[4];
  const ScalarType one_div_delta_wnd = parameters_real_[5];

  t_idx = (t - t0) * one_div_delta;
  const int t_wnd_idx = (get_tau(t_idx) - t - t0_window) * one_div_delta_wnd;
  delta_t = get_tau(t_idx) - t - get_window_tau(t_wnd_idx);

  const int oversampling = parameters_int_[3];
  const int window_sampling = parameters_int_[4];
  t_idx -= oversampling;
  const int delta_tau_index = t_wnd_idx - oversampling * window_sampling;
  const int j = delta_tau_index % window_sampling;
  const int i = (delta_tau_index - j) / window_sampling;

  constexpr int interpolation_size = mode == CUBIC ? 4 : 2;
  conv_coeff_idx = interpolation_size * i + interpolation_size * 4 * oversampling * j;
}

template <typename ScalarType>
__host__ NfftHelperManager<ScalarType>::~NfftHelperManager() {
  cudaFree(NfftHelper<ScalarType>::sub_matrix_);
  cudaFree(NfftHelper<ScalarType>::parameters_int_);
  cudaFree(NfftHelper<ScalarType>::parameters_real_);
}

}  // details
}  // nfft
}  // math
}  // dca

#endif  // DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER
