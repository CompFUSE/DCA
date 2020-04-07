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
// Integration tests for the Dnfft1D class.

#include "dca/math/nfft/dnfft_1d.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
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
void computeWithDnfft(
    const std::vector<typename DnfftType::ElementType>& t,
    const std::vector<typename DnfftType::ElementType>& f, DnfftType& dnfft_obj,
    function<std::complex<typename DnfftType::ElementType>, dmn_variadic<FreqDmn, OtherDmn>>& f_w);

template <typename Real>
void computeWithDft(const std::vector<Real>& t, const std::vector<Real>& f,
                    function<std::complex<Real>, dmn_variadic<FreqDmn, OtherDmn>>& f_w);

void initializeDomains(double beta, int time_slices) {
  static bool initialized = false;

  if (!initialized) {
    dca::phys::domains::time_domain::initialize(beta, time_slices, 1.e-16);
    const int positive_frequencies = time_slices - 2;
    dca::phys::domains::frequency_domain::initialize(beta, positive_frequencies);
    initialized = true;
  }
}

template <typename Real>
class Dnfft1DTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(Dnfft1DTest, TestTypes);

TYPED_TEST(Dnfft1DTest, CubicInterpolation) {
  using Real = TypeParam;

  // Initialize time and frequency domains.
  const Real beta = 10.;
  const int time_slices = 100;
  initializeDomains(beta, time_slices);

  // Prepare random samples.
  dca::math::random::StdRandomWrapper<std::mt19937_64> rng(0, 1, 0);
  const int samples = 1e4;

  std::vector<Real> t(samples);
  std::vector<Real> f(samples);

  const double begin = TimeDmn::get_elements().front();
  const double delta = TimeDmn::get_elements().back() - TimeDmn::get_elements().front();

  for (int l = 0; l < samples; ++l) {
    const Real t_val = begin + rng() * delta;
    t[l] = t_val;
    f[l] = std::exp(-2. * M_PI / delta * (t_val - begin));
  }

  // Compute f(w) using the discrete Fourier transform (DFT).
  function<std::complex<Real>, dmn_variadic<FreqDmn, OtherDmn>> f_w_dft("f_w_dft");
  computeWithDft(t, f, f_w_dft);

  // Compute f(w) using the delayed-NFFT algorithm.
  constexpr int oversampling = 8;
  dca::math::nfft::Dnfft1D<Real, FreqDmn, OtherDmn, oversampling, dca::math::nfft::CUBIC> dnfft_obj;

  function<std::complex<Real>, dmn_variadic<FreqDmn, OtherDmn>> f_w_dnfft("f_w_dnfft");
  computeWithDnfft(t, f, dnfft_obj, f_w_dnfft);

  // Check errors.
  const auto err = dca::func::util::difference(f_w_dft, f_w_dnfft);

  const Real tolerance = std::is_same<Real, double>::value ? 1e-9 : 1e-5;
  EXPECT_LT(err.l1, tolerance);
  EXPECT_LT(err.l2, tolerance);
  EXPECT_LT(err.l_inf, tolerance);
}

template <typename DnfftType>
void computeWithDnfft(
    const std::vector<typename DnfftType::ElementType>& t,
    const std::vector<typename DnfftType::ElementType>& f, DnfftType& dnfft_obj,
    function<std::complex<typename DnfftType::ElementType>, dmn_variadic<FreqDmn, OtherDmn>>& f_w) {
  dnfft_obj.resetAccumulation();

  const double begin = TimeDmn::get_elements().front();
  const double delta = TimeDmn::get_elements().back() - TimeDmn::get_elements().front();

  for (int t_ind = 0; t_ind < t.size(); ++t_ind) {
    // Transform the time point to the interval [-0.5, 0.5].
    const double scaled_t = (t[t_ind] - begin) / delta - 0.5;
    dnfft_obj.accumulate(0, scaled_t, f[t_ind]);
  }

  dnfft_obj.finalize(f_w);
}

template <typename Real>
void computeWithDft(const std::vector<Real>& t, const std::vector<Real>& f,
                    function<std::complex<Real>, dmn_variadic<FreqDmn, OtherDmn>>& f_w) {
  const std::complex<Real> i(0, 1);

  f_w = 0.;

  for (int o_ind = 0; o_ind < OtherDmn::dmn_size(); ++o_ind) {
    for (int t_ind = 0; t_ind < t.size(); ++t_ind) {
      for (int w_ind = 0; w_ind < FreqDmn::dmn_size(); ++w_ind) {
        const Real t_val = t[t_ind];
        const Real w_val = FreqDmn::get_elements()[w_ind];

        f_w(w_ind, o_ind) += f[t_ind] * std::exp(i * t_val * w_val);
      }
    }
  }
}
