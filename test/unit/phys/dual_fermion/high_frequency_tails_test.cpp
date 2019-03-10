// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests high_frequency_tails.hpp

#include "dca/phys/dual_fermion/high_frequency_tails.hpp"

#include <complex>
#include <limits>

#include "gtest/gtest.h"

#include "dca/parallel/no_concurrency/no_concurrency.hpp"

using namespace dca;

namespace testing {
// testing::

// Mock parameters class for initialization of the frequency domains.
struct MockParameters {
  double get_beta() const {
    return 10.;
  }
  int get_sp_fermionic_frequencies() const {
    return 128;
  }
  int get_four_point_fermionic_frequencies() const {
    return 64;
  }
};
}  // namespace testing

class HighFrequencyTailsTest : public ::testing::Test {
protected:
  using SpFreqDmn = phys::df::HighFrequencyTails::SpFreqDmn;
  using TpFreqDmn = phys::df::HighFrequencyTails::TpFreqDmn;

  using OtherDmns =
      func::dmn_variadic<func::dmn_0<func::dmn<3, int>>, func::dmn_0<func::dmn<2, double>>>;

  HighFrequencyTailsTest() : concurrency_(0, nullptr) {}

  static void SetUpTestCase() {
    SpFreqDmn::parameter_type::initialize(parameters_);
    TpFreqDmn::parameter_type::initialize(parameters_);
  }

  static const testing::MockParameters parameters_;

  const parallel::NoConcurrency concurrency_;

  func::function<std::complex<double>, func::dmn_variadic<OtherDmns, SpFreqDmn>> Sigma_sp_freq_;
  func::function<std::complex<double>, func::dmn_variadic<OtherDmns, TpFreqDmn>> Sigma_tp_freq_;
};

const testing::MockParameters HighFrequencyTailsTest::parameters_;

TEST_F(HighFrequencyTailsTest, ExactFit) {
  const double A = 2.8;
  const double B = 3.2;

  auto fitting_func = [A, B](const double w, const int shift) {
    return std::complex<double>((B + shift) / (w * w), -(A + shift) / w);
  };

  // Prepare input.
  for (int w_ind = 0; w_ind < TpFreqDmn::dmn_size(); ++w_ind) {
    const auto w = TpFreqDmn::get_elements()[w_ind];
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      Sigma_tp_freq_(o_ind, w_ind) = fitting_func(w, o_ind);
    }
  }

  const int tail_freqs = 32;

  phys::df::HighFrequencyTails::compute(concurrency_, tail_freqs, Sigma_tp_freq_, Sigma_sp_freq_);

  const auto tol = std::numeric_limits<double>::epsilon();

  for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
    const auto w = SpFreqDmn::get_elements()[w_ind];
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      const std::complex<double> expected = fitting_func(w, o_ind);
      EXPECT_NEAR(expected.real(), Sigma_sp_freq_(o_ind, w_ind).real(), tol);
      EXPECT_NEAR(expected.imag(), Sigma_sp_freq_(o_ind, w_ind).imag(), tol);
    }
  }
}

TEST_F(HighFrequencyTailsTest, InvalidTailFreqs) {
  const int tail_freqs = TpFreqDmn::dmn_size() / 2 + 1;

  EXPECT_THROW(phys::df::HighFrequencyTails::compute(concurrency_, tail_freqs, Sigma_tp_freq_,
                                                     Sigma_sp_freq_),
               std::invalid_argument);
}
