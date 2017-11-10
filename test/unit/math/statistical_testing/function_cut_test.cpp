// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the utility provided by function_cut.hpp.

#include "dca/math/statistical_testing/function_cut.hpp"

#include "gtest/gtest.h"

std::string test_folder = DCA_SOURCE_DIR "/test/unit/math/statistical_testing/";

TEST(FunctionCutTest, FrequencyCut) {
  // Initialize the physics domains.
  const int n_b = 2, n_k = 4, n_frq = 8;
  // Initialize band domain.
  dca::math::util::details::Bdmn::parameter_type::get_size() = n_b;
  // Initialize mock momentum space domain.
  using Kdmn = dca::func::dmn_0<dca::func::dmn<n_k, double>>;
  // Initialize frequency domain
  dca::math::util::details::Wdmn::parameter_type::initialize(1., n_frq);

  // Write a mock Sigma function.
  dca::func::function<std::complex<double>, dca::math::util::SigmaDomain<Kdmn>> sigma;
  const int n_spin_band_sqr = 4 * n_b * n_b;
  EXPECT_EQ(n_spin_band_sqr * n_k * 2 * n_frq, sigma.size());
  for (int k = 0; k < n_k; ++k) {
    int frq_val = -n_frq;
    for (int w = 0; w < 2 * n_frq; ++w) {
      sigma(0, 0, 0, 0, k, w) = std::complex<double>(k, frq_val);
      ++frq_val;
      if (frq_val == 0)
        ++frq_val;
    }
  }

  const int n_frq_kept = 2;
  // Compute cut function: only n_frq_kept positive frequencies, the diagonal of the up sector, and
  // the K points are kept. Real and imaginary parts are unrolled into two real entries.
  auto sigma_cut = dca::math::util::cutFrequency(sigma, n_frq_kept);
  EXPECT_EQ(n_b * n_b * n_k * 2 * n_frq_kept, sigma_cut.size());

  for (int k = 0; k < n_k; ++k)
    for (int w = 0; w < n_frq_kept; ++w)
      for (int re_im = 0; re_im < 2; ++re_im) {
        const double expected = re_im == 0 ? k : w + 1;
        EXPECT_DOUBLE_EQ(expected, sigma_cut(0, 0, k, w, re_im));
      }
}

TEST(FunctionCutTest, BandDiagonal) {
  using dca::func::dmn;
  using dca::func::dmn_0;
  using dca::func::dmn_variadic;
  using MockWDmn = dmn_0<dmn<8, double>>;
  using MockKDmn = dmn_0<dmn<4, double>>;
  using Bdmn = dca::math::util::details::Bdmn;
  Bdmn::parameter_type::get_size() = 3;

  dca::func::function<int, dmn_variadic<Bdmn, Bdmn, MockKDmn, MockWDmn>> f_initial;
  for (int i = 0; i < f_initial.size(); ++i)
    f_initial(i) = i;

  dca::func::function<int, dmn_variadic<Bdmn, MockKDmn, MockWDmn>> f_diagonal =
      dca::math::util::bandDiagonal(f_initial);

  for (int b = 0; b < Bdmn::dmn_size(); ++b)
    for (int k = 0; k < MockKDmn::dmn_size(); ++k)
      for (int w = 0; w < MockWDmn::dmn_size(); ++w)
        EXPECT_EQ(f_initial(b, b, k, w), f_diagonal(b, k, w));
}
