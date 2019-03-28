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
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"

using namespace dca;

class BareDualGreensFunctionTest : public ::testing::Test {
protected:
  using Complex = std::complex<double>;

  using SpinDmn = func::dmn_0<func::dmn<2, int>>;
  using BandDmn = func::dmn_0<func::dmn<1, int>>;
  using BandSpinDmn = func::dmn_variadic<BandDmn, SpinDmn>;
  using KClusterDmn = func::dmn_0<func::dmn<2, std::vector<double>>>;
  using KSuperlatticeDmn = func::dmn_0<func::dmn<4, std::vector<double>>>;
  using MatsubaraFreqDmn = func::dmn_0<func::dmn<2, double>>;

  using BareDualGreensFunctionType =
      phys::df::BareDualGreensFunction<Complex, parallel::NoConcurrency, BandSpinDmn, KClusterDmn,
                                       KSuperlatticeDmn, MatsubaraFreqDmn>;

  BareDualGreensFunctionTest()
      : concurrency_(0, nullptr), bare_dual_gf_(concurrency_, eps_bar_, eps_hat_, Delta_, G_) {}

  const parallel::NoConcurrency concurrency_;

  BareDualGreensFunctionType::Dispersion eps_bar_;
  BareDualGreensFunctionType::NonTransInvariantDispersion eps_hat_;
  BareDualGreensFunctionType::GF Delta_;
  BareDualGreensFunctionType::GF G_;

  BareDualGreensFunctionType bare_dual_gf_;
};

TEST_F(BareDualGreensFunctionTest, TwoSitesDiagonalInOmegaKtilde) {
  // Prepare input.
  const Complex eps_bar_0(4.2, 2.4);
  const Complex eps_bar_1(3.7, 7.3);
  for (int s = 0; s < BandSpinDmn::dmn_size(); ++s) {
    eps_bar_(s, s, 0) = eps_bar_0;
    eps_bar_(s, s, 1) = eps_bar_1;
  }

  const Complex eps_hat_00(1., 1.);
  const Complex eps_hat_01(2., 4.);
  const Complex eps_hat_10(3., 9.);
  const Complex eps_hat_11(4., 16.);
  for (int k = 0; k < KSuperlatticeDmn::dmn_size(); ++k) {
    eps_hat_(0, 0, k) = eps_hat_00;
    eps_hat_(0, 1, k) = eps_hat_01;
    eps_hat_(1, 0, k) = eps_hat_10;
    eps_hat_(1, 1, k) = eps_hat_11;
  }

  const Complex Delta_0(1.23, -1.23);
  const Complex Delta_1(4.56, 4.56);
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < BandSpinDmn::dmn_size(); ++s) {
      Delta_(s, s, 0, w) = Delta_0;
      Delta_(s, s, 1, w) = Delta_1;
    }

  const Complex G_0(3.14, 1.59);
  const Complex G_1(2.65, 3.59);
  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w)
    for (int s = 0; s < BandSpinDmn::dmn_size(); ++s) {
      G_(s, s, 0, w) = G_0;
      G_(s, s, 1, w) = G_1;
    }

  // Compute expected result.
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

  bare_dual_gf_.compute();
  const auto& G0_tilde = bare_dual_gf_.get();

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
