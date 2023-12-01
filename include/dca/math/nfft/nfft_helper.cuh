// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// Helper class for computing interpolation indices during the single particle accumulation.

#ifndef DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER_CUH
#define DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER_CUH

#include <array>
#include <cassert>
#include "dca/platform/dca_gpu.h"
#include <stdexcept>

#include "dca/config/mc_options.hpp"
#include "dca/math/nfft/nfft_mode_names.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/cluster_helper.cuh"

namespace dca {
namespace math {
namespace nfft {
namespace details {
// dca::math::nfft::details::

template <typename REAL>
class NfftHelper {
public:
  using Real = REAL;

  static void set(int nb, const int nc, const int* add_r, int lda, const int* sub_r, int lds, Real t0,
                  Real delta_t, const Real t0_window, const Real delta_t_window, const Real beta);

  __device__ int inline computeLinearIndex(int b1, int b2, int r1, int r2) const;

  __device__ auto inline computeTau(Real t1, Real t2) const;

  __device__ inline int deltaR(int r1, int r2) const;

  template <NfftModeNames mode, int oversampling, int window_sampling, typename TType>
  __device__ inline void computeInterpolationIndices(TType t, int& t_idx, int& conv_coeff_idx,
                                                     TType& delta_t) const;

private:
  int nb_;
  Real t0_, delta_t_, t0_window_, delta_t_window_, one_div_delta_, one_div_delta_t_window_,
      tau_scale_;
};

/** Global instances
 *  type tagged
 */
extern __device__ __constant__ NfftHelper<double> nfft_helper_double;
extern __device__ __constant__ NfftHelper<float> nfft_helper_float;
extern template void NfftHelper<double>::set(int nb, const int nc, const int* add_r, int lda,
                                             const int* sub_r, int lds, double t0, double delta_t,
                                             const double t0_window, const double delta_t_window,
                                             const double beta);
extern template void NfftHelper<float>::set(int nb, const int nc, const int* add_r, int lda,
                                            const int* sub_r, int lds, float t0, float delta_t,
                                            const float t0_window, const float delta_t_window,
                                            const float beta);

template <typename REAL>
__device__ int inline NfftHelper<REAL>::computeLinearIndex(const int b1, const int b2, const int r1,
                                                           const int r2) const {
  const int delta_r = dca::phys::solver::details::cluster_real_helper.subtract(r2, r1);
  return b1 + nb_ * (b2 + nb_ * delta_r);
}

template <typename REAL>
__device__ auto NfftHelper<REAL>::computeTau(REAL t1, REAL t2) const {
  return (t1 - t2) * tau_scale_;
}

template <typename REAL>
template <NfftModeNames mode, int oversampling, int window_sampling, typename TType>
__device__ void NfftHelper<REAL>::computeInterpolationIndices(TType t, int& t_idx,
                                                              int& conv_coeff_idx,
                                                              TType& delta_t) const {
  static_assert(mode == CUBIC || mode == LINEAR,
                "This method is defined only for linear or cubic interpolations.");

  auto get_tau = [=](int index) { return t0_ + index * delta_t_; };
  auto get_window_tau = [&](int index) { return t0_window_ + index * delta_t_window_; };

  t_idx = (t - t0_) * one_div_delta_;
  const int t_wnd_idx = (get_tau(t_idx) - t - t0_window_) * one_div_delta_t_window_;
  delta_t = get_tau(t_idx) - t - get_window_tau(t_wnd_idx);

  t_idx -= oversampling;
  const int delta_tau_index = t_wnd_idx - oversampling * window_sampling;
  const int i = delta_tau_index / window_sampling;
  const int j = delta_tau_index - i * window_sampling;

  constexpr int interpolation_size = mode == CUBIC ? 4 : 2;
  conv_coeff_idx = interpolation_size * i + interpolation_size * 4 * oversampling * j;
}

}  // namespace details
}  // namespace nfft
}  // namespace math
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_MATH_NFFT_NFFT_HELPER_CUH
