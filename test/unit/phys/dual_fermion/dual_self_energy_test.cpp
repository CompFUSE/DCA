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
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"

using namespace dca;

class DualSelfEnergyTest : public ::testing::Test {
protected:
  static constexpr int dimension_ = 2;
  static constexpr int num_pos_tp_freqs_ = 3;
  static constexpr int num_exchange_freqs_ = 2;

  using BandDmn = func::dmn_0<func::dmn<1, int>>;

  using TpFreqDmn = func::dmn_0<func::dmn<2 * num_pos_tp_freqs_, double>>;
  using FreqExchangeDmn = func::dmn_0<func::dmn<num_exchange_freqs_, int>>;

  using RClusterDmn = phys::ClusterDomainAliases<dimension_>::RClusterDmn;
  using KClusterDmn = phys::ClusterDomainAliases<dimension_>::KClusterDmn;

  using RSuperlatticeDmn = phys::ClusterDomainAliases<dimension_>::RSpSuperlatticeDmn;
  using KSuperlatticeDmn = phys::ClusterDomainAliases<dimension_>::KSpSuperlatticeDmn;

  using DualSelfEnergyType =
      phys::df::DualSelfEnergy<double, parallel::NoConcurrency, BandDmn, KClusterDmn,
                               KSuperlatticeDmn, TpFreqDmn, FreqExchangeDmn>;

  using TpGreensFunction = DualSelfEnergyType::TpGreensFunction;
  using DualGreensFunction = DualSelfEnergyType::DualGreensFunction;

  DualSelfEnergyTest()
      : concurrency_(0, nullptr),
        beta_(10.),
        Nc_(RClusterDmn::dmn_size()),
        V_(RSuperlatticeDmn::dmn_size()),
        Gamma_long_uu_val_(3.14),
        Gamma_long_ud_val_(42.),
        Gamma_trans_ud_val_(1.23),
        Sigma_tilde_(concurrency_, beta_, G0_tilde_, Gamma_long_uu_, Gamma_long_ud_, Gamma_trans_ud_) {
  }

  static void SetUpTestCase() {
    const std::array<double, 4> cluster_basis{1., 0., 0., 1.};
    const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 1}};
    phys::domains::cluster_domain_initializer<RClusterDmn>::execute(cluster_basis.data(),
                                                                    cluster_superbasis);

    const std::array<double, 4> superlattice_basis{2., 0., 0., 2.};
    const int superlattice_size = 4;
    const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                                {0, superlattice_size}};
    phys::domains::cluster_domain_initializer<RSuperlatticeDmn>::execute(superlattice_basis.data(),
                                                                         superlattice_superbasis);
  }

  const parallel::NoConcurrency concurrency_;

  const double beta_;

  const int Nc_;
  const int V_;

  DualGreensFunction G0_tilde_;

  const double Gamma_long_uu_val_;
  TpGreensFunction Gamma_long_uu_;

  const double Gamma_long_ud_val_;
  TpGreensFunction Gamma_long_ud_;

  const double Gamma_trans_ud_val_;
  TpGreensFunction Gamma_trans_ud_;

  DualSelfEnergyType Sigma_tilde_;
};

TEST_F(DualSelfEnergyTest, Compute1stOrder) {
  // Prepare G0_tilde_ w/o w or k_tilde depedency.
  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      G0_tilde_(0, 0, k_tilde, w) = 1.;
      G0_tilde_(0, 1, k_tilde, w) = 2.;
      G0_tilde_(0, 1, k_tilde, w) = 3.;
      G0_tilde_(1, 1, k_tilde, w) = 4.;
    }

  // Prepare Gamma_long_uu_ and Gamma_long_ud_ w/o K1, K2, K_ex, w1 or w2 dependecy.
  for (int w_ex = 0; w_ex < FreqExchangeDmn::dmn_size(); ++w_ex)
    for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2)
      for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1)
        for (int K_ex = 0; K_ex < KClusterDmn::dmn_size(); ++K_ex)
          for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
            for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
              // Choose different value for w_ex == 0 because it is the only frequency used.
              if (w_ex == 0) {
                Gamma_long_uu_(0, 0, 0, 0, K1, K2, K_ex, w1, w2, w_ex) = Gamma_long_uu_val_;
                Gamma_long_ud_(0, 0, 0, 0, K1, K2, K_ex, w1, w2, w_ex) = Gamma_long_ud_val_;
              }
              else {
                Gamma_long_uu_(0, 0, 0, 0, K1, K2, K_ex, w1, w2, w_ex) = -1.;
                Gamma_long_ud_(0, 0, 0, 0, K1, K2, K_ex, w1, w2, w_ex) = -1.;
              }
            }

  Sigma_tilde_.compute1stOrder();
  const DualGreensFunction& Sigma_tilde_1st = Sigma_tilde_.get();

  const double prefactor = -1. / (Nc_ * V_ * beta_) * (Gamma_long_uu_val_ + Gamma_long_ud_val_) *
                           KSuperlatticeDmn::dmn_size() * TpFreqDmn::dmn_size();

  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      EXPECT_DOUBLE_EQ(
          (prefactor * (G0_tilde_(0, 0, k_tilde, w) + G0_tilde_(1, 1, k_tilde, w))).real(),
          Sigma_tilde_1st(0, 0, k_tilde, w).real());

      EXPECT_DOUBLE_EQ(
          (prefactor * (G0_tilde_(0, 1, k_tilde, w) + G0_tilde_(1, 0, k_tilde, w))).real(),
          Sigma_tilde_1st(0, 1, k_tilde, w).real());

      EXPECT_DOUBLE_EQ(
          (prefactor * (G0_tilde_(1, 0, k_tilde, w) + G0_tilde_(0, 1, k_tilde, w))).real(),
          Sigma_tilde_1st(1, 0, k_tilde, w).real());

      EXPECT_DOUBLE_EQ(
          (prefactor * (G0_tilde_(1, 1, k_tilde, w) + G0_tilde_(0, 0, k_tilde, w))).real(),
          Sigma_tilde_1st(1, 1, k_tilde, w).real());
    }

  // Imaginary part should be zero since all input functions are real.
  for (int w = 0; w < TpFreqDmn::dmn_size(); ++w)
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde)
      for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2)
        for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1)
          EXPECT_EQ(0., Sigma_tilde_1st(K1, K2, k_tilde, w).imag());
}

TEST_F(DualSelfEnergyTest, Compute2ndOrder) {
  Sigma_tilde_.compute2ndOrder();
  const DualGreensFunction& Sigma_tilde_2nd = Sigma_tilde_.get();
}
