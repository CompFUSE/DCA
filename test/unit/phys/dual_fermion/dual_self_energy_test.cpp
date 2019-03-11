// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dual_self_energy.hpp.

#include "dca/phys/dual_fermion/dual_self_energy.hpp"

#include <array>
#include <complex>
#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

using namespace dca;

namespace testing {
// testing::

// Mock parameters class for initialization of the frequency domains.
struct MockParameters {
  double get_beta() const {
    return 10.;
  }
  int get_sp_fermionic_frequencies() const {
    return 8;
  }
  int get_four_point_fermionic_frequencies() const {
    return 2;
  }
  bool compute_all_transfers() const {
    return true;
  }
  int get_four_point_frequency_transfer() const {
    return 2;
  }
};
}  // namespace testing

class DualSelfEnergyTest : public ::testing::Test {
protected:
  static constexpr int dimension_ = 2;

  using DualSelfEnergyType = phys::df::DualSelfEnergy<double, parallel::NoConcurrency, dimension_>;

  using BandDmn = DualSelfEnergyType::BandDmn;

  using SpFreqDmn = func::dmn_0<phys::domains::frequency_domain>;
  using FreqExchangeDmn = DualSelfEnergyType::FreqExchangeDmn;
  using TpFreqDmn = DualSelfEnergyType::TpFreqDmn;
  using ExtFreqDmn = DualSelfEnergyType::ExtFreqDmn;

  using KClusterDmn = DualSelfEnergyType::KClusterDmn;
  using KSuperlatticeDmn = DualSelfEnergyType::KSuperlatticeDmn;

  using TpGF = DualSelfEnergyType::TpGF;
  using DualGFExtFreq = DualSelfEnergyType::DualGFExtFreq;
  using DualGFTpFreq = DualSelfEnergyType::DualGFTpFreq;

  DualSelfEnergyTest()
      : concurrency_(0, nullptr),
        Gamma_long_uu_val_(3.14),
        Gamma_long_ud_val_(42.),
        Gamma_tran_ud_val_(1.23),
        Sigma_tilde_(concurrency_, parameters_.get_beta(), G0_tilde_, Gamma_long_uu_,
                     Gamma_long_ud_, Gamma_tran_ud_) {}

  static void SetUpTestCase() {
    const int num_bands = 1;
    BandDmn::parameter_type::initialize(parameters_, num_bands, std::vector<int>(num_bands),
                                        std::vector<std::vector<double>>(num_bands));

    SpFreqDmn::parameter_type::initialize(parameters_);
    FreqExchangeDmn::parameter_type::initialize(parameters_);
    TpFreqDmn::parameter_type::initialize(parameters_);
    ExtFreqDmn::parameter_type::initialize(parameters_);

    const std::array<double, 4> cluster_basis{1., 0., 0., 1.};
    const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 1}};
    phys::domains::cluster_domain_initializer<func::dmn_0<KClusterDmn::parameter_type::dual_type>>::execute(
        cluster_basis.data(), cluster_superbasis);

    const std::array<double, 4> superlattice_basis{2., 0., 0., 1.};
    const int superlattice_size = 2;
    const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                                {0, superlattice_size}};
    phys::domains::cluster_domain_initializer<func::dmn_0<KSuperlatticeDmn::parameter_type::dual_type>>::execute(
        superlattice_basis.data(), superlattice_superbasis);
  }

  static const testing::MockParameters parameters_;

  const parallel::NoConcurrency concurrency_;

  DualGFExtFreq G0_tilde_;

  const double Gamma_long_uu_val_;
  TpGF Gamma_long_uu_;

  const double Gamma_long_ud_val_;
  TpGF Gamma_long_ud_;

  const double Gamma_tran_ud_val_;
  TpGF Gamma_tran_ud_;

  DualSelfEnergyType Sigma_tilde_;
};

const testing::MockParameters DualSelfEnergyTest::parameters_;

TEST_F(DualSelfEnergyTest, Compute1stOrder) {
  // Prepare G0_tilde_ w/o w or k_tilde depedency.
  for (int w = 0; w < ExtFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      G0_tilde_(0, 0, k_tilde, w) = 1.;
      G0_tilde_(0, 1, k_tilde, w) = 2.;
      G0_tilde_(1, 0, k_tilde, w) = 3.;
      G0_tilde_(1, 1, k_tilde, w) = 4.;
    }

  // Prepare Gamma_long_uu_ and Gamma_long_ud_ w/o K1, K2, Q, wn or wm dependecy.
  for (int l = 0; l < FreqExchangeDmn::dmn_size(); ++l)
    for (int m = 0; m < TpFreqDmn::dmn_size(); ++m)
      for (int n = 0; n < TpFreqDmn::dmn_size(); ++n)
        for (int Q = 0; Q < KClusterDmn::dmn_size(); ++Q)
          for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
            for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
              // Choose different value for l == 0 because it is the only frequency used.
              if (l == 0) {
                Gamma_long_uu_(0, 0, 0, 0, K1, K2, Q, n, m, l) = Gamma_long_uu_val_;
                Gamma_long_ud_(0, 0, 0, 0, K1, K2, Q, n, m, l) = Gamma_long_ud_val_;
              }
              else {
                Gamma_long_uu_(0, 0, 0, 0, K1, K2, Q, n, m, l) = -1.;
                Gamma_long_ud_(0, 0, 0, 0, K1, K2, Q, n, m, l) = -1.;
              }
            }

  Sigma_tilde_.compute1stOrder();
  const DualGFTpFreq& Sigma_tilde_1st = Sigma_tilde_.get();

  const double prefactor =
      -1. / (KClusterDmn::dmn_size() * KSuperlatticeDmn::dmn_size() * parameters_.get_beta()) *
      (Gamma_long_uu_val_ + Gamma_long_ud_val_) * KSuperlatticeDmn::dmn_size() *
      TpFreqDmn::dmn_size();

  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      EXPECT_NEAR((prefactor * (G0_tilde_(0, 0, 0, 0) + G0_tilde_(1, 1, 0, 0))).real(),
                  Sigma_tilde_1st(0, 0, k_tilde, w).real(),
                  1000 * std::numeric_limits<double>::epsilon());

      EXPECT_NEAR((prefactor * (G0_tilde_(0, 1, 0, 0) + G0_tilde_(1, 0, 0, 0))).real(),
                  Sigma_tilde_1st(0, 1, k_tilde, w).real(),
                  1000 * std::numeric_limits<double>::epsilon());

      EXPECT_NEAR((prefactor * (G0_tilde_(1, 0, 0, 0) + G0_tilde_(0, 1, 0, 0))).real(),
                  Sigma_tilde_1st(1, 0, k_tilde, w).real(),
                  1000 * std::numeric_limits<double>::epsilon());

      EXPECT_NEAR((prefactor * (G0_tilde_(1, 1, 0, 0) + G0_tilde_(0, 0, 0, 0))).real(),
                  Sigma_tilde_1st(1, 1, k_tilde, w).real(),
                  1000 * std::numeric_limits<double>::epsilon());
    }

  // Imaginary part should be zero since all input functions are real.
  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde)
      for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
        for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1)
          EXPECT_EQ(0., Sigma_tilde_1st(K1, K2, k_tilde, w).imag());
}

TEST_F(DualSelfEnergyTest, Compute2ndOrder) {
  // Prepare G0_tilde_ w/o w or k_tilde depedency.
  for (int w = 0; w < ExtFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      G0_tilde_(0, 0, k_tilde, w) = 1.;
      G0_tilde_(0, 1, k_tilde, w) = 2.;
      G0_tilde_(1, 0, k_tilde, w) = 3.;
      G0_tilde_(1, 1, k_tilde, w) = 4.;
    }

  // Prepare constant Gamma's.
  Gamma_long_uu_ = Gamma_long_uu_val_;
  Gamma_long_ud_ = Gamma_long_ud_val_;
  Gamma_tran_ud_ = Gamma_tran_ud_val_;

  Sigma_tilde_.compute2ndOrder();
  const DualGFTpFreq& Sigma_tilde_2nd = Sigma_tilde_.get();

  const double prefactor =
      -1. /
      (2. * KClusterDmn::dmn_size() * KClusterDmn::dmn_size() * KSuperlatticeDmn::dmn_size() *
       KSuperlatticeDmn::dmn_size() * parameters_.get_beta() * parameters_.get_beta()) *
      (Gamma_long_uu_val_ * Gamma_long_uu_val_ + Gamma_long_ud_val_ * Gamma_long_ud_val_ +
       Gamma_tran_ud_val_ * Gamma_tran_ud_val_) *
      KSuperlatticeDmn::dmn_size() * KSuperlatticeDmn::dmn_size() * TpFreqDmn::dmn_size() *
      (2 * FreqExchangeDmn::dmn_size() - 1);

  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      // K1 = 0, K2 = 0
      EXPECT_NEAR(
          (prefactor * (G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0)))
              .real(),
          Sigma_tilde_2nd(0, 0, k_tilde, w).real(),
          10000000 * std::numeric_limits<double>::epsilon());

      // K1 = 0, K2 = 1
      EXPECT_NEAR(
          (prefactor * (G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0)))
              .real(),
          Sigma_tilde_2nd(0, 1, k_tilde, w).real(),
          10000000 * std::numeric_limits<double>::epsilon());

      // K1 = 1, K2 = 0
      EXPECT_NEAR(
          (prefactor * (G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0)))
              .real(),
          Sigma_tilde_2nd(1, 0, k_tilde, w).real(),
          10000000 * std::numeric_limits<double>::epsilon());

      // K1 = 1, K2 = 1
      EXPECT_NEAR(
          (prefactor * (G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 1, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(1, 0, 0, 0) +
                        G0_tilde_(0, 0, 0, 0) * G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(0, 1, 0, 0) * G0_tilde_(1, 0, 0, 0) * G0_tilde_(0, 0, 0, 0) +
                        G0_tilde_(1, 1, 0, 0) * G0_tilde_(0, 0, 0, 0) * G0_tilde_(0, 0, 0, 0)))
              .real(),
          Sigma_tilde_2nd(1, 1, k_tilde, w).real(),
          10000000 * std::numeric_limits<double>::epsilon());
    }

  // Imaginary part should be zero since all input functions are real.
  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde)
      for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
        for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1)
          EXPECT_EQ(0., Sigma_tilde_2nd(K1, K2, k_tilde, w).imag());
}
