// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class implements the 1D delayed-NFFT (d-NFFT) algorithm [1].
// It requires an FFTW library with the FFTW3 interface.
//
// References:
// [1] P. Staar, T. A. Maier, and T. C. Schulthess, J. Phys.: Conf. Ser. 402, 012015 (2012).

#ifndef DCA_MATH_NFFT_DNFFT_1D_HPP
#define DCA_MATH_NFFT_DNFFT_1D_HPP

#include <algorithm>
#include <cassert>
#include <complex>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <fftw3.h>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/nfft/domains/domains.hpp"
#include "dca/math/nfft/nfft_atomic_convolution.hpp"
#include "dca/math/nfft/nfft_mode_names.hpp"
#include "dca/math/nfft/window_functions/gaussian_window_function.hpp"
#include "dca/math/nfft/window_functions/kaiser_bessel_function.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

// The method to compute the kernel is determined by mode.
// Options are EXACT [evaluation], LINEAR [interpolation], and CUBIC [interpolation] (recommended).
template <typename ScalarType, typename WDmn, typename PDmn, int oversampling = 8, NfftModeNames mode = CUBIC>
class Dnfft1D {
public:
  using ThisType = Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>;
  using ElementType = ScalarType;

  Dnfft1D();
  Dnfft1D(ThisType&& other) = default;

  void resetAccumulation();

  // Adds the sample (t_val, f_val) to the accumulated function.
  // linind is the linear index of the sample w.r.t p_dmn (= all non-transformed (discrete)
  // domains).
  // Preconditions: t_val must be in the interval [-0.5, 0.5].
  void accumulate(int linind, ScalarType t_val, ScalarType f_val);
  // Version with subindices.
  // subind contains the subindices of the sample w.r.t. the subdomains of p_dmn.
  void accumulate(const int* subind, ScalarType t_val, ScalarType f_val);

  // Performs the final FFT on the accumulated function.
  // Out: f_w
  template <typename OtherScalarType>
  void finalize(func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w);

  // Sums the accumulated data in the time domain.
  ThisType& operator+=(const ThisType& other_one);

  constexpr int get_oversampling() const {
    return oversampling;
  }
  constexpr int get_window_sampling() const {
    return window_sampling_;
  }

  int maximumFrequency() const {
    return WDmn::dmn_size() / 2;
  }

protected:
  static constexpr int window_sampling_ = 32;
  static constexpr double window_function_sigma_ = 2.;

  using WindowFunction = kaiser_bessel_function;

  using LinearCoefficientsDmn = func::dmn_0<nfft_linear_coefficients_domain>;
  using CubicCoefficientsDmn = func::dmn_0<nfft_cubic_coefficients_domain>;

  using OversamplingDmn = func::dmn_0<nfft_oversampling_domain<ThisType>>;
  using WindowSamplingDmn = func::dmn_0<nfft_window_sampling_domain<ThisType>>;

  using PaddedTimeDmn = func::dmn_0<nfft_time_domain<PADDED, ThisType>>;
  using LeftOrientedTimeDmn = func::dmn_0<nfft_time_domain<LEFT_ORIENTED, ThisType>>;
  using WindowFunctionTimeDmn = func::dmn_0<nfft_time_domain<WINDOW_FUNCTION, ThisType>>;

  using ConvolutionTimeDmn = func::dmn_variadic<OversamplingDmn, WindowSamplingDmn>;

  using PaddedTimePDmn = func::dmn_variadic<PaddedTimeDmn, PDmn>;
  using LeftOrientedPDmn = func::dmn_variadic<LeftOrientedTimeDmn, PDmn>;

  static inline auto& get_convolution_time_values();
  static inline auto& get_linear_convolution_matrices();
  static inline auto& get_cubic_convolution_matrices();

  func::function<ScalarType, PaddedTimePDmn> f_tau_;

private:
  static void initializeDomains(const ThisType& this_obj);
  static void initializeStaticFunctions();

  void convoluteToFTauExact(int index, ScalarType t_val, ScalarType f_val);
  void convoluteToFTauFineLinearInterpolation(int index, ScalarType t_val, ScalarType f_val);
  void convoluteToFTauFineCubicInterpolation(int index, ScalarType t_val, ScalarType f_val);

  void foldTimeDomainBack();

  template <typename OtherScalarType>
  void transformFTauToFW(
      func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w) const;

  static inline ScalarType tau(int idx) {
    return PaddedTimeDmn::get_elements()[idx];
  }
  static inline ScalarType fineTau(int idx) {
    return WindowFunctionTimeDmn::get_elements()[idx];
  }

  static inline auto& get_phi_wn();
};

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::Dnfft1D() : f_tau_("f_tau_") {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    initializeDomains(*this);
    initializeStaticFunctions();
  });
  f_tau_.reset();
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::resetAccumulation() {
  f_tau_ = 0.;
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::initializeDomains(const ThisType& this_obj) {
  OversamplingDmn::parameter_type::initialize(this_obj);
  WindowSamplingDmn::parameter_type::initialize(this_obj);

  nfft_time_domain<LEFT_ORIENTED, ThisType>::initialize(this_obj);
  nfft_time_domain<PADDED, ThisType>::initialize(this_obj);
  nfft_time_domain<WINDOW_FUNCTION, ThisType>::initialize(this_obj);
  nfft_time_domain<FOLDED_WINDOW_FUNCTION, ThisType>::initialize(this_obj);

  if (mode == LINEAR || mode == CUBIC) {
    auto& convolution_times = get_convolution_time_values();
    convolution_times.reset();

    for (int i = 0; i < OversamplingDmn::dmn_size(); ++i)
      for (int j = 0; j < WindowSamplingDmn::dmn_size(); ++j)
        convolution_times(i, j) =
            nfft_time_domain<WINDOW_FUNCTION,
                             ThisType>::get_elements()[j + i * WindowSamplingDmn::dmn_size()];
  }
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::initializeStaticFunctions() {
  WindowFunction::n = PaddedTimeDmn::dmn_size();
  WindowFunction::m = oversampling;

  WindowFunction::sigma = window_function_sigma_;

  if (mode == LINEAR)
    get_linear_convolution_matrices().reset();
  else if (mode == CUBIC)
    get_cubic_convolution_matrices().reset();

  if (mode == LINEAR || mode == CUBIC) {
    int index = 0;
    auto& convolution_time_values = get_convolution_time_values();
    ScalarType delta = convolution_time_values(0, 1) - convolution_time_values(0, 0);

    for (int i = 0; i < OversamplingDmn::dmn_size(); ++i) {
      for (int j = 0; j < WindowSamplingDmn::dmn_size(); ++j) {
        assert(std::abs(convolution_time_values(i, j) - fineTau(index)) < 1.e-6);

        ScalarType tau = convolution_time_values(i, j);

        ScalarType f0 = WindowFunction::phi_t(tau);
        ScalarType f1 = WindowFunction::phi_t(tau + delta);

        ScalarType df0 = WindowFunction::d_phi_t(tau);
        ScalarType df1 = WindowFunction::d_phi_t(tau + delta);

        ScalarType a = f0;
        ScalarType b = df0;

        ScalarType c = -(3. * f0 - 3. * f1 + 2. * df0 * delta + df1 * delta) / std::pow(delta, 2);
        ScalarType d = -(-2. * f0 + 2. * f1 - 1. * df0 * delta - df1 * delta) / std::pow(delta, 3);

        if (mode == LINEAR) {
          get_linear_convolution_matrices()(0, i, j) = f0;
          get_linear_convolution_matrices()(1, i, j) = (f1 - f0) / delta;
        }
        else if (mode == CUBIC) {
          get_cubic_convolution_matrices()(0, i, j) = a;
          get_cubic_convolution_matrices()(1, i, j) = b;
          get_cubic_convolution_matrices()(2, i, j) = c;
          get_cubic_convolution_matrices()(3, i, j) = d;
        }
        index += 1;
      }
    }
  }

  auto& phi_wn = get_phi_wn();
  const auto& matsubara_freq_indices = WDmn::parameter_type::get_indices();
  for (int l = 0; l < WDmn::dmn_size(); ++l)
    phi_wn(l) = WindowFunction::phi_wn(matsubara_freq_indices[l]);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
inline void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::accumulate(const int linind,
                                                                            const ScalarType t_val,
                                                                            const ScalarType f_val) {
  assert(t_val > -0.5 - 1.e-6 && t_val < 0.5 + 1.e-6);

  switch (mode) {
    case EXACT:
      convoluteToFTauExact(linind, t_val, f_val);
      break;

    case LINEAR:
      convoluteToFTauFineLinearInterpolation(linind, t_val, f_val);
      break;

    case CUBIC:
      convoluteToFTauFineCubicInterpolation(linind, t_val, f_val);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
inline void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::accumulate(const int* const subind,
                                                                            const ScalarType t_val,
                                                                            const ScalarType f_val) {
  const static PDmn p_dmn_obj;
  std::size_t linind = 0;
  p_dmn_obj.subind_2_linind(subind, linind);
  accumulate(linind, t_val, f_val);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
template <typename OtherScalarType>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::finalize(
    func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w) {
  foldTimeDomainBack();
  transformFTauToFW(f_w);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::convoluteToFTauExact(
    const int index, const ScalarType t_val, const ScalarType f_val) {
  assert(t_val > -0.5 - 1.e-6 && t_val < 0.5 + 1.e-6);

  const ScalarType T_0 = PaddedTimeDmn::parameter_type::first_element();
  const ScalarType one_div_Delta = PaddedTimeDmn::parameter_type::get_one_div_Delta();

  int lambda_0 = (t_val - T_0) * one_div_Delta;

  for (int l = -oversampling; l <= oversampling; ++l)
    f_tau_(lambda_0 + l, index) += f_val * WindowFunction::phi_t(tau(lambda_0 + l) - t_val);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
inline void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::convoluteToFTauFineLinearInterpolation(
    const int index, const ScalarType t_val, const ScalarType f_val) {
  assert(t_val > -0.5 - 1.e-6 && t_val < 0.5 + 1.e-6);

  const ScalarType t_0 = WindowFunctionTimeDmn::parameter_type::first_element();
  const ScalarType T_0 = PaddedTimeDmn::parameter_type::first_element();

  const ScalarType one_div_delta = PaddedTimeDmn::parameter_type::get_one_div_delta();
  const ScalarType one_div_Delta = PaddedTimeDmn::parameter_type::get_one_div_Delta();

  int tau_0 = (t_val - T_0) * one_div_Delta;
  int tau_1 = (tau(tau_0) - t_val - t_0) * one_div_delta;

  assert(tau(tau_0) - 1.e-10 < t_val && t_val < tau(tau_0 + 1) + 1.e-10);

  ScalarType diff_tau = tau(tau_0) - t_val - fineTau(tau_1);

  assert(diff_tau > -1.e-6 && diff_tau < PaddedTimeDmn::parameter_type::get_delta());

  ScalarType y_ptr[2];

  y_ptr[0] = f_val;
  y_ptr[1] = f_val * diff_tau;

  int tau_index = tau_0 - oversampling;
  int delta_tau_index = tau_1 - oversampling * window_sampling_;

  int j = delta_tau_index % window_sampling_;
  int i = (delta_tau_index - j) / window_sampling_;
  assert(delta_tau_index == i * window_sampling_ + j);

  ScalarType* f_tau_ptr = &f_tau_(tau_index, index);
  ScalarType* matrix_ptr = &get_linear_convolution_matrices()(0, i, j);

  nfft_atomic_convolution<2 * oversampling + 1, 0>::execute_linear(f_tau_ptr, matrix_ptr, y_ptr);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
inline void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::convoluteToFTauFineCubicInterpolation(
    const int index, const ScalarType t_val, const ScalarType f_val) {
  assert(t_val > -0.5 - 1.e-6 && t_val < 0.5 + 1.e-6);

  const ScalarType t_0 = WindowFunctionTimeDmn::parameter_type::first_element();
  const ScalarType T_0 = PaddedTimeDmn::parameter_type::first_element();

  const ScalarType delta = PaddedTimeDmn::parameter_type::get_delta();
  const ScalarType Delta = PaddedTimeDmn::parameter_type::get_Delta();

  const ScalarType one_div_delta = PaddedTimeDmn::parameter_type::get_one_div_delta();
  const ScalarType one_div_Delta = PaddedTimeDmn::parameter_type::get_one_div_Delta();

  int tau_0 = (t_val - T_0) * one_div_Delta;
  ScalarType t0_val_lb = T_0 + tau_0 * Delta;

  assert(tau(tau_0) - 1.e-6 < t_val && t_val < tau(tau_0 + 1) + 1.e-6);
  assert(std::abs(tau(tau_0) - t0_val_lb) < 1.e-6);

  // int tau_1 = (tau(tau_0)-t_val-t_0)*one_div_delta;
  int tau_1 = (t0_val_lb - t_val - t_0) * one_div_delta;
  ScalarType t1_val_lb = t_0 + tau_1 * delta;

  ScalarType diff_tau = t0_val_lb - t_val - t1_val_lb;  // fineTau(tau_1);

  assert(diff_tau > -1.e-6 && diff_tau < PaddedTimeDmn::parameter_type::get_delta());

  ScalarType y_ptr[4];

  y_ptr[0] = f_val;
  y_ptr[1] = y_ptr[0] * diff_tau;
  y_ptr[2] = y_ptr[1] * diff_tau;
  y_ptr[3] = y_ptr[2] * diff_tau;

  int tau_index = tau_0 - oversampling;
  int delta_tau_index = tau_1 - oversampling * window_sampling_;

  int j = delta_tau_index % window_sampling_;
  int i = (delta_tau_index - j) / window_sampling_;
  assert(delta_tau_index == i * window_sampling_ + j);

  ScalarType* f_tau_ptr = &f_tau_(tau_index, index);
  ScalarType* matrix_ptr = &get_cubic_convolution_matrices()(0, i, j);

  nfft_atomic_convolution<2 * oversampling + 1, 0>::execute_cubic(f_tau_ptr, matrix_ptr, y_ptr);
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::foldTimeDomainBack() {
  // Fold the halos of f_tau_ to enforce periodicity of the time domain.

  const int n_padded = nfft_time_domain<PADDED, ThisType>::get_size();
  const int n = nfft_time_domain<LEFT_ORIENTED, ThisType>::get_size();
  const int padding = (n_padded - n) / 2;
  assert(n >= padding);

  for (int p_ind = 0; p_ind < PDmn::dmn_size(); ++p_ind) {
    // Fold left side.
    for (int t_ind = 0; t_ind < padding; ++t_ind)
      f_tau_(t_ind + n, p_ind) += f_tau_(t_ind, p_ind);
    // Fold right side.
    for (int t_ind = n_padded - padding; t_ind < n_padded; ++t_ind)
      f_tau_(t_ind - n, p_ind) += f_tau_(t_ind, p_ind);
  }
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
template <typename OtherScalarType>
void Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::transformFTauToFW(
    func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w) const {
  const int n_padded = nfft_time_domain<PADDED, ThisType>::get_size();
  const int n = nfft_time_domain<LEFT_ORIENTED, ThisType>::get_size();
  const int padding = (n_padded - n) / 2;

  std::vector<double> f_in(n);

  // We need to use std::complex<double>, because fftw_complex is just a typedef for double[2] and
  // as such does not meet the requirements of "Erasable" required by std::vector.
  std::vector<std::complex<double>> f_out(n / 2 + 1);

  static std::mutex fftw_mutex;
  fftw_mutex.lock();
  // See http://www.fftw.org/fftw3_doc/Complex-numbers.html for why the cast should be safe.
  fftw_plan plan = fftw_plan_dft_r2c_1d(
      n, f_in.data(), reinterpret_cast<fftw_complex*>(f_out.data()), FFTW_ESTIMATE);
  fftw_mutex.unlock();

  func::function<std::complex<ScalarType>, LeftOrientedPDmn> f_omega;
  for (int p_ind = 0; p_ind < PDmn::dmn_size(); ++p_ind) {
    std::copy_n(&f_tau_(padding, p_ind), n, f_in.data());

    fftw_execute(plan);

    for (int t_ind = 0; t_ind < n / 2; ++t_ind) {
      f_omega(t_ind, p_ind).real(-f_out[t_ind].real());
      f_omega(t_ind, p_ind).imag(f_out[t_ind].imag());
    }

    for (int t_ind = n / 2; t_ind < n; ++t_ind) {
      f_omega(t_ind, p_ind).real(-f_out[n - t_ind].real());
      f_omega(t_ind, p_ind).imag(-f_out[n - t_ind].imag());
    }
  }

  fftw_destroy_plan(plan);

  std::vector<int> w_indices(0);
  const auto& matsubara_freq_indices = WDmn::parameter_type::get_indices();
  for (int w_ind = 0; w_ind < WDmn::dmn_size(); ++w_ind) {
    for (int t_ind = 0; t_ind < n; ++t_ind) {
      if (matsubara_freq_indices[w_ind] == t_ind || matsubara_freq_indices[w_ind] + n == t_ind) {
        w_indices.push_back(t_ind);
        break;
      }
    }
  }

  for (int p_ind = 0; p_ind < PDmn::dmn_size(); ++p_ind)
    for (int w_ind = 0; w_ind < WDmn::dmn_size(); ++w_ind)
      f_w(w_ind, p_ind) = f_omega(w_indices[w_ind], p_ind) / get_phi_wn()(w_ind);

  f_w *= 1. / n;
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>& Dnfft1D<ScalarType, WDmn, PDmn, oversampling,
                                                             mode>::operator+=(const ThisType& other) {
  f_tau_ += other.f_tau_;
  return *this;
}

// Static member getters.
template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
auto& Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::get_convolution_time_values() {
  static func::function<ScalarType, func::dmn_variadic<OversamplingDmn, WindowSamplingDmn>>
      convolution_time_values("convolution_times");
  return convolution_time_values;
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
auto& Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::get_linear_convolution_matrices() {
  static func::function<ScalarType,
                        func::dmn_variadic<LinearCoefficientsDmn, OversamplingDmn, WindowSamplingDmn>>
      linear_convolution_matrices("linear_convolution_matrices");
  return linear_convolution_matrices;
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
auto& Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::get_cubic_convolution_matrices() {
  static func::function<ScalarType,
                        func::dmn_variadic<CubicCoefficientsDmn, OversamplingDmn, WindowSamplingDmn>>
      cubic_convolution_matrices("cubic_convolution_matrices");
  return cubic_convolution_matrices;
}

template <typename ScalarType, typename WDmn, typename PDmn, int oversampling, NfftModeNames mode>
auto& Dnfft1D<ScalarType, WDmn, PDmn, oversampling, mode>::get_phi_wn() {
  static func::function<ScalarType, WDmn> phi_wn("phi_wn");
  return phi_wn;
}

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_DNFFT_1D_HPP
