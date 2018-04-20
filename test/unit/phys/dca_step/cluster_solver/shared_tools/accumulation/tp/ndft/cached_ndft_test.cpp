// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// Integration tests for the cached_ndft class.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/profiling/events/time.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/dnft_test.hpp"

constexpr int n_samples = 40;
constexpr int n_bands = 2;
constexpr int n_frqs = 2;
using TestSetup = dca::testing::DnftTest<n_samples, n_bands, n_frqs>;

double computeWithFastDNFT(const TestSetup::Configuration& config, const TestSetup::Matrix& M,
                           TestSetup::F_w_w& f_w);
void computeWithDft(const TestSetup::Configuration& config, const TestSetup::Matrix& M,
                    TestSetup::F_w_w& f_w);

TEST_F(TestSetup, Execute) {
  // Compute the DNFT with the CachedNdft class.
  F_w_w f_w_fast("f_w_fast");
  const double time = computeWithFastDNFT(configuration_, M_, f_w_fast);

  const auto err = dca::func::util::difference(f_baseline_, f_w_fast);
  EXPECT_LT(err.l_inf, 1e-14);

  std::cout << "\nCached ndft time [sec]:\t " << time << "\n";
}

double computeWithFastDNFT(const TestSetup::Configuration& config, const TestSetup::Matrix& M,
                           TestSetup::F_w_w& f_w) {
  dca::func::function<std::complex<double>,
                      dca::func::dmn_variadic<TestSetup::BDmn, TestSetup::BDmn, TestSetup::RDmn,
                                              TestSetup::RDmn, TestSetup::FreqPosDmn, TestSetup::FreqDmn>>
      f_b_b_r_r_w_w;
  dca::phys::solver::accumulator::CachedNdft<double, TestSetup::RDmn, TestSetup::FreqDmn,
                                             TestSetup::FreqPosDmn, dca::linalg::CPU>
      nft_obj;

  dca::profiling::WallTime start_time;
  nft_obj.execute(config, M, f_b_b_r_r_w_w);
  dca::profiling::WallTime end_time;

  // Rearrange output.
  const int n_w = TestSetup::FreqPosDmn::dmn_size();
  auto invert_w = [=](const int w) { return 2 * n_w - 1 - w; };
  for (int b2 = 0; b2 < TestSetup::BDmn::dmn_size(); ++b2)
    for (int b1 = 0; b1 < TestSetup::BDmn::dmn_size(); ++b1)
      for (int r2 = 0; r2 < TestSetup::RDmn::dmn_size(); ++r2)
        for (int r1 = 0; r1 < TestSetup::RDmn::dmn_size(); ++r1)
          for (int w2 = 0; w2 < TestSetup::FreqDmn::dmn_size(); ++w2)
            for (int w1 = 0; w1 < n_w; ++w1) {
              f_w(b1, b2, r1, r2, w1 + n_w, w2) = f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2);
              f_w(b1, b2, r1, r2, invert_w(w1 + n_w), invert_w(w2)) =
                  std::conj(f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2));
            }

  dca::profiling::Duration duration(end_time, start_time);
  return duration.sec + 1.e-6 * duration.usec;
}
