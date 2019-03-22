// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests bare_dual_greens_function.hpp.

#include "dca/phys/dual_fermion/bare_dual_greens_function.hpp"

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
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn>> eps_bar;
  const std::complex<double> eps_bar_0(4.2, 2.4);
  const std::complex<double> eps_bar_1(3.7, 7.3);
  for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
    eps_bar(s, s, 0) = eps_bar_0;
    eps_bar(s, s, 1) = eps_bar_1;
  }

  func::function<std::complex<double>, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn>> eps_hat;
  const std::complex<double> eps_hat_00(1., 1.);
  const std::complex<double> eps_hat_01(2., 4.);
  const std::complex<double> eps_hat_10(3., 9.);
  const std::complex<double> eps_hat_11(4., 16.);
  for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k) {
    eps_hat(0, 0, k) = eps_hat_00;
    eps_hat(0, 1, k) = eps_hat_01;
    eps_hat(1, 0, k) = eps_hat_10;
    eps_hat(1, 1, k) = eps_hat_11;
  }

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>
      Delta;
  const std::complex<double> Delta_0(1.23, -1.23);
  const std::complex<double> Delta_1(4.56, 4.56);
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
      Delta(s, s, 0, w) = Delta_0;
      Delta(s, s, 1, w) = Delta_1;
    }

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>
      G;
  const std::complex<double> G_0(3.14, 1.59);
  const std::complex<double> G_1(2.65, 3.59);
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < OrbitalSpinDmn::dmn_size(); ++s) {
      G(s, s, 0, w) = G_0;
      G(s, s, 1, w) = G_1;
    }

  func::function<std::complex<double>,
                 func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, MatsubaraFreqDmn>>
      G0_tilde;

  phys::df::computeBareDualGreensFunction(eps_bar, eps_hat, Delta, G, concurrency_, G0_tilde);

  const auto det_D_inv =
      1. / ((Delta_0 + eps_bar_0 - eps_hat_00) * (Delta_1 + eps_bar_1 - eps_hat_11) -
            eps_hat_01 * eps_hat_10);

  const auto det_A_inv = 1. / ((G_0 + det_D_inv * (Delta_1 + eps_bar_1 - eps_hat_11)) *
                                   (G_1 + det_D_inv * (Delta_0 + eps_bar_0 - eps_hat_00)) -
                               det_D_inv * det_D_inv * eps_hat_01 * eps_hat_10);

  const auto G0_tilde_00 =
      -det_A_inv * G_0 * G_0 * (G_1 + det_D_inv * (Delta_0 + eps_bar_0 - eps_hat_00));

  const auto G0_tilde_11 =
      -det_A_inv * G_1 * G_1 * (G_0 + det_D_inv * (Delta_1 + eps_bar_1 - eps_hat_11));

  const auto G0_tilde_01 = det_A_inv * G_0 * G_1 * det_D_inv * eps_hat_01;

  const auto G0_tilde_10 = det_A_inv * G_1 * G_0 * det_D_inv * eps_hat_10;

  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_tilde_00.real(), G0_tilde(0, 0, k, w).real());
      EXPECT_DOUBLE_EQ(G0_tilde_00.imag(), G0_tilde(0, 0, k, w).imag());

      EXPECT_DOUBLE_EQ(G0_tilde_01.real(), G0_tilde(0, 1, k, w).real());
      EXPECT_DOUBLE_EQ(G0_tilde_01.imag(), G0_tilde(0, 1, k, w).imag());

      EXPECT_DOUBLE_EQ(G0_tilde_10.real(), G0_tilde(1, 0, k, w).real());
      EXPECT_DOUBLE_EQ(G0_tilde_10.imag(), G0_tilde(1, 0, k, w).imag());

      EXPECT_DOUBLE_EQ(G0_tilde_11.real(), G0_tilde(1, 1, k, w).real());
      EXPECT_DOUBLE_EQ(G0_tilde_11.imag(), G0_tilde(1, 1, k, w).imag());
    }
}
