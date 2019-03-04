// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_bare_dual_greens_function.hpp.

#include "dca/phys/dual_fermion/compute_bare_dual_greens_function.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"

using namespace dca;

class ComputeBareDualGreensFunctionTest : public ::testing::Test {
protected:
  using SpinDmn = func::dmn_0<func::dmn<2, int>>;
  using OrbitalDmn = func::dmn_0<func::dmn<1, int>>;
  using OrbitalSpinDmn = func::dmn_variadic<OrbitalDmn, SpinDmn>;
  using KClusterDmn = func::dmn_0<func::dmn<2, std::vector<double>>>;
  using KSuperlatticeDmn = func::dmn_0<func::dmn<4, std::vector<double>>>;
  using MatsubaraFreqDmn = func::dmn_0<func::dmn<2, double>>;

  ComputeBareDualGreensFunctionTest() : concurrency_(0, nullptr) {}

  const parallel::NoConcurrency concurrency_;
};

TEST_F(ComputeBareDualGreensFunctionTest, TwoClusterSites) {
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn>> eps_cg;
  const double eps_cg_0 = 42.;
  const double eps_cg_1 = 99.;
  for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
    eps_cg(s, s, 0) = eps_cg_0;
    eps_cg(s, s, 1) = eps_cg_1;
  }

  func::function<std::complex<double>, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn>>
      eps_tilde;
  for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k) {
    eps_tilde(0, 0, k) = 1.;
    eps_tilde(0, 1, k) = 2.;
    eps_tilde(1, 0, k) = 3.;
    eps_tilde(1, 1, k) = 4.;
  }

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>
      Delta;
  const double Delta_0 = 1.23;
  const double Delta_1 = 4.56;
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
      Delta(s, s, 0, w) = Delta_0;
      Delta(s, s, 1, w) = Delta_1;
    }

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>
      G;
  const double G_0 = 3.14;
  const double G_1 = 2.72;
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
      G(s, s, 0, w) = G_0;
      G(s, s, 1, w) = G_1;
    }

  func::function<std::complex<double>,
                 func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, MatsubaraFreqDmn>>

      G0_tilde;

  phys::df::computeBareDualGreensFunction(eps_cg, eps_tilde, Delta, G, concurrency_, G0_tilde);

  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k) {
      EXPECT_DOUBLE_EQ(
          -G_0 * 1. / (G_0 + 1. / ((Delta_0 + eps_cg_0) - eps_tilde(0, 0, k).real())) * G_0,
          G0_tilde(0, 0, k, w).real());
      EXPECT_DOUBLE_EQ(-G_0 * 1. / (G_1 + 1. / (-eps_tilde(0, 1, k).real())) * G_1,
                       G0_tilde(0, 1, k, w).real());
      EXPECT_DOUBLE_EQ(-G_1 * 1. / (G_0 + 1. / (-eps_tilde(1, 0, k).real())) * G_0,
                       G0_tilde(1, 0, k, w).real());
      EXPECT_DOUBLE_EQ(
          -G_1 * 1. / (G_1 + 1. / ((Delta_1 + eps_cg_1) - eps_tilde(1, 1, k).real())) * G_1,
          G0_tilde(1, 1, k, w).real());
    }

  // Imaginary part should be zero since all input functions are real.
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1)
        for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k)
          EXPECT_EQ(0., G0_tilde(K1, K2, k, w).imag());
}
