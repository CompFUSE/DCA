// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests high_frequency_tails.hpp for the case when SpFreqDmn and TpFreqDmn have the same
// size.

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
    return get_sp_fermionic_frequencies();
  }
};
}  // namespace testing

class HighFrequencyTailsTest : public ::testing::Test {
protected:
  using OtherDmns =
      func::dmn_variadic<func::dmn_0<func::dmn<3, int>>, func::dmn_0<func::dmn<2, double>>>;

  using HighFrequencyTailsType =
      phys::df::HighFrequencyTails<double, parallel::NoConcurrency, OtherDmns>;

  using SpFreqDmn = HighFrequencyTailsType::SpFreqDmn;
  using TpFreqDmn = HighFrequencyTailsType::TpFreqDmn;

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
  // Prepare input.
  for (int w_ind = 0; w_ind < TpFreqDmn::dmn_size(); ++w_ind) {
    const auto w = TpFreqDmn::get_elements()[w_ind];
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      Sigma_tp_freq_(o_ind, w_ind) = w + o_ind;
    }
  }

  const int tail_freqs = 32;
  HighFrequencyTailsType high_freq_tails(concurrency_, tail_freqs, Sigma_tp_freq_, Sigma_sp_freq_);

  high_freq_tails.compute();

  for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      EXPECT_DOUBLE_EQ(Sigma_tp_freq_(o_ind, w_ind).real(), Sigma_sp_freq_(o_ind, w_ind).real());
      EXPECT_DOUBLE_EQ(Sigma_tp_freq_(o_ind, w_ind).imag(), Sigma_sp_freq_(o_ind, w_ind).imag());
    }
  }
}
