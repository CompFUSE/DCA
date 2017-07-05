// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Integration tests for the Dnfft1D class.

#include "dca/math/nfft/dnfft_1d.hpp"

#include <algorithm>  // std::max
#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

using dca::func::function;
using dca::func::dmn_variadic;

using TimeDmn = dca::func::dmn_0<dca::phys::domains::time_domain>;
using FreqDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;

// Represents all non-transformed domains.
using OtherDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

template <typename DnfftType>
void computeWithDnfft(const std::vector<double>& t, const std::vector<double>& f,
                      DnfftType& dnfft_obj,
                      function<std::complex<double>, dmn_variadic<FreqDmn, OtherDmn>>& f_w);
void computeWithDft(const std::vector<double>& t, const std::vector<double>& f,
                    function<std::complex<double>, FreqDmn>& f_w);

TEST(Dnfft1DTest, CubicInterpolation) {
  // Initialize time and frequency domains.
  const double beta = 10.;
  const int time_slices = 100;
  const int positive_frequencies = time_slices - 2;

  dca::phys::domains::time_domain::initialize(beta, time_slices, 1.e-16);
  dca::phys::domains::frequency_domain::initialize(beta, positive_frequencies);

  // Prepare random samples.
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);
  const int samples = 1e4;

  std::vector<double> t(samples);
  std::vector<double> f(samples);

  const double begin = TimeDmn::get_elements().front();
  const double delta = TimeDmn::get_elements().back() - TimeDmn::get_elements().front();

  for (int l = 0; l < samples; ++l) {
    const double t_val = begin + rng() * delta;
    t[l] = t_val;
    f[l] = std::exp(-2. * M_PI / delta * (t_val - begin));
  }

  // Compute f(w) using the discrete Fourier transform (DFT).
  function<std::complex<double>, FreqDmn> f_w_dft("f_w_dft");
  computeWithDft(t, f, f_w_dft);

  // Compute f(w) using the delayed-NFFT algorithm.
  constexpr int oversampling = 8;
  dca::math::nfft::Dnfft1D<double, FreqDmn, OtherDmn, oversampling, dca::math::nfft::CUBIC> dnfft_obj;

  function<std::complex<double>, dmn_variadic<FreqDmn, OtherDmn>> f_w_dnfft("f_w_dnfft");
  computeWithDnfft(t, f, dnfft_obj, f_w_dnfft);

  // Check errors.
  double f_w_dft_l1 = 0.;
  double f_w_dft_l2 = 0.;
  double f_w_dft_linf = 0.;

  double l1_error = 0.;
  double l2_error = 0.;
  double linf_error = 0.;

  for (int w_ind = 0; w_ind < FreqDmn::dmn_size(); w_ind++) {
    const double f_w_dft_abs = std::abs(f_w_dft(w_ind));
    f_w_dft_l1 += f_w_dft_abs;
    f_w_dft_l2 += f_w_dft_abs * f_w_dft_abs;
    f_w_dft_linf = std::max(f_w_dft_linf, f_w_dft_abs);

    const double err = std::abs(f_w_dft(w_ind) - f_w_dnfft(w_ind, 0));
    l1_error += err;
    l2_error += err * err;
    linf_error = std::max(linf_error, err);
  }

  l1_error /= f_w_dft_l1;
  l2_error = std::sqrt(l2_error / f_w_dft_l2);
  linf_error /= f_w_dft_linf;

  // std::cout << "l1_error = " << l1_error << std::endl;
  // std::cout << "l2_error = " << l2_error << std::endl;
  // std::cout << "linf_error = " << linf_error << std::endl;

  EXPECT_LT(l1_error, 1.e-9);
  EXPECT_LT(l2_error, 1.e-9);
  EXPECT_LT(linf_error, 1.e-9);
}

template <typename DnfftType>
void computeWithDnfft(const std::vector<double>& t, const std::vector<double>& f,
                      DnfftType& dnfft_obj,
                      function<std::complex<double>, dmn_variadic<FreqDmn, OtherDmn>>& f_w) {
  dnfft_obj.initialize();

  const double begin = TimeDmn::get_elements().front();
  const double delta = TimeDmn::get_elements().back() - TimeDmn::get_elements().front();

  for (int t_ind = 0; t_ind < t.size(); ++t_ind) {
    // Transform the time point to the interval [-0.5, 0.5].
    const double scaled_t = (t[t_ind] - begin) / delta - 0.5;
    dnfft_obj.accumulate(0, scaled_t, f[t_ind]);
  }

  dnfft_obj.finalize(f_w);
}

void computeWithDft(const std::vector<double>& t, const std::vector<double>& f,
                    function<std::complex<double>, FreqDmn>& f_w) {
  const std::complex<double> i(0, 1);

  f_w = 0.;

  for (int t_ind = 0; t_ind < t.size(); ++t_ind) {
    for (int w_ind = 0; w_ind < FreqDmn::dmn_size(); ++w_ind) {
      const double t_val = t[t_ind];
      const double w_val = FreqDmn::get_elements()[w_ind];

      f_w(w_ind) += f[t_ind] * std::exp(i * t_val * w_val);
    }
  }
}
