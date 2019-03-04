// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_hybridization_function.hpp.

#include "dca/phys/dca_algorithms/compute_hybridization_function.hpp"

#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"

using namespace dca;

class ComputeHybridizationFunctionTest : public ::testing::Test {
protected:
  using SpinDmn = func::dmn<2, int>;
  using OrbitalDmn = func::dmn<1, int>;  // 1 orbital
  using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<OrbitalDmn>, func::dmn_0<SpinDmn>>;
  using KClusterDmn = func::dmn<1, std::vector<double>>;
  using MatsubaraFreqDmn = func::dmn<4, double>;  // Matsubara frequency domain with 4 elements.

  ComputeHybridizationFunctionTest() : concurrency_(0, nullptr) {}

  static void SetUpTestCase() {
    const std::vector<std::vector<double>> k_vecs{{0., 0.}};
    KClusterDmn::set_elements(k_vecs);

    const double beta = 1.;
    const std::vector<double> freqs{-3. * M_PI / beta, -M_PI / beta, M_PI / beta, 3. * M_PI / beta};
    MatsubaraFreqDmn::set_elements(freqs);
  }

  const parallel::NoConcurrency concurrency_;
};

TEST_F(ComputeHybridizationFunctionTest, SingleSiteSquareLattice) {
  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>>>
      eps_K_cg;
  // Set diagonal elements (in orbital-spin domain) of eps_K_cg to some constant.
  const double eps_K_cg_val = -5.2;
  eps_K_cg(0, 0, 0, 0, 0) = eps_K_cg(0, 1, 0, 1, 0) = eps_K_cg_val;

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>,
                                    func::dmn_0<MatsubaraFreqDmn>>>
      G0_cluster_K_wn;

  // Set diagonal elements (in orbital-spin domain) of G0_cluster_K_wn to some constant.
  const double G0_cluster_K_wn_val = 8.4;
  for (int wn = 0; wn < MatsubaraFreqDmn::dmn_size(); ++wn)
    G0_cluster_K_wn(0, 0, 0, 0, 0, wn) = G0_cluster_K_wn(0, 1, 0, 1, 0, wn) = G0_cluster_K_wn_val;

  const double mu = 0.5;  // chemical potential

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KClusterDmn>,
                                    func::dmn_0<MatsubaraFreqDmn>>>
      Delta_K_wn;

  phys::computeHybridizationFunction(eps_K_cg, G0_cluster_K_wn, mu, concurrency_, Delta_K_wn);

  for (int wn = 0; wn < MatsubaraFreqDmn::dmn_size(); ++wn) {
    // Off-diagonal elements in orbital-spin domain should be zero.
    EXPECT_DOUBLE_EQ(0., Delta_K_wn(0, 0, 0, 1, 0, wn).real());
    EXPECT_DOUBLE_EQ(0., Delta_K_wn(0, 0, 0, 1, 0, wn).imag());
    EXPECT_DOUBLE_EQ(0., Delta_K_wn(0, 1, 0, 0, 0, wn).real());
    EXPECT_DOUBLE_EQ(0., Delta_K_wn(0, 1, 0, 0, 0, wn).imag());

    // Check diagonal elements.
    const double expected_real = mu - eps_K_cg_val - 1. / G0_cluster_K_wn_val;
    EXPECT_DOUBLE_EQ(expected_real, Delta_K_wn(0, 0, 0, 0, 0, wn).real());
    EXPECT_DOUBLE_EQ(expected_real, Delta_K_wn(0, 1, 0, 1, 0, wn).real());

    const double expected_imag = MatsubaraFreqDmn::get_elements()[wn];
    EXPECT_DOUBLE_EQ(expected_imag, Delta_K_wn(0, 0, 0, 0, 0, wn).imag());
    EXPECT_DOUBLE_EQ(expected_imag, Delta_K_wn(0, 1, 0, 1, 0, wn).imag());
  }
}
