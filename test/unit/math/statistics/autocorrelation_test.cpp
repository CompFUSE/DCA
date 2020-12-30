// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests autocorrelation.hpp by generating a time series of known autocorrelation time.

#include "dca/math/statistics/autocorrelation.hpp"

#include <random>
#include "gtest/gtest.h"

TEST(AutocorrelationTest, All) {
  const std::vector<double> tau_exp_vals{1, 5, 3};
  const std::vector<double> means{0, 1, -1};
  const std::vector<double> stdevs{1, 0.5, 2};

  const int n_samples = 5e4;

  std::mt19937_64 rng(42);

  for (int test_id = 0; test_id < tau_exp_vals.size(); ++test_id) {
    const double tau_exp = tau_exp_vals[test_id];

    std::vector<double> independent_series(n_samples);
    std::normal_distribution<double> distro(means[test_id], stdevs[test_id]);

    const int n_corr = 100;
    dca::math::statistics::Autocorrelation<double> autocorr(n_corr);

    // Compute and accumulate time series.
    for (auto& x : independent_series) {
      x = distro(rng);
    }

    for (int i = n_corr; i < independent_series.size(); ++i) {
      double correlated = 0;

      for (int j = i - n_corr; j <= i; ++j) {
        correlated += std::exp(-(i - j) / tau_exp) * independent_series[j];
      }

      autocorr.addSample(correlated);
    }

    // See: http://www.hep.fsu.edu/~berg/teach/mcmc08/material/lecture07mcmc3.pdf
    const auto tau_integrated = 1. + (2 * std::exp(-1. / tau_exp)) / (1 - std::exp(-1. / tau_exp));

    const double tau_computed = autocorr.computeAutocorrelationTime();

    const double rel_diff = 2 * (tau_computed - tau_integrated) / (tau_computed + tau_integrated);
    EXPECT_NEAR(rel_diff, 0, 0.1);
  }
}

TEST(AutocorrelationTest, ComplexCompiles) {
  dca::math::statistics::Autocorrelation<std::complex<float>> corr(10);

  corr.addSample(std::complex<float>(1, 0));
  corr.addSample(std::complex<float>(1, -1));

  // Note: this computation will fail to converge.
  corr.computeAutocorrelationTime();
}
