// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_reducible_tp_vertex.hpp.

#include "dca/phys/dca_algorithms/compute_reducible_tp_vertex.hpp"

#include <array>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/four_point_type.hpp"

using namespace dca;

class ComputeReducibleTpVertexTest : public ::testing::Test {
protected:
  static constexpr int dimension_ = 2;
  static constexpr double beta_ = 10.;
  static constexpr int Nc_ = 4;
  static constexpr int num_pos_tp_freqs_ = 3;
  static constexpr int num_exchange_freqs_ = 2;

  using BandDmn = func::dmn_0<func::dmn<1, int>>;
  using SpinDmn = func::dmn_0<func::dmn<2, int>>;

  using SpFreqDmn = func::dmn_0<func::dmn<2 * (num_pos_tp_freqs_ + num_exchange_freqs_), double>>;
  using TpFreqDmn = func::dmn_0<func::dmn<2 * num_pos_tp_freqs_, double>>;
  using FreqExchangeDmn = func::dmn_0<func::dmn<num_exchange_freqs_, int>>;

  using RDmn = phys::ClusterDomainAliases<dimension_>::RClusterDmn;
  using KDmn = phys::ClusterDomainAliases<dimension_>::KClusterDmn;

  using SpDomain = func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, KDmn, SpFreqDmn>;
  using SpGreensFunction = func::function<std::complex<double>, SpDomain>;

  using TpDomain = func::dmn_variadic<BandDmn, BandDmn, BandDmn, BandDmn, KDmn, KDmn, KDmn,
                                      TpFreqDmn, TpFreqDmn, FreqExchangeDmn>;
  using TpGreensFunction = func::function<std::complex<double>, TpDomain>;

  ComputeReducibleTpVertexTest()
      : concurrency_(0, nullptr),
        G_val_(3.14, 99),
        G4_val_(42, -2.72),
        G4_over_GGGG(beta_ * Nc_ * G4_val_ / (G_val_ * G_val_ * G_val_ * G_val_)),
        one_over_GG(beta_ * Nc_ / (G_val_ * G_val_)) {
    for (int w = 0; w < SpFreqDmn::dmn_size(); ++w)
      for (int k = 0; k < KDmn::dmn_size(); ++k)
        G_(0, 0, 0, 0, k, w) = G_val_;

    for (int w_ex = 0; w_ex < FreqExchangeDmn::dmn_size(); ++w_ex)
      for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2)
        for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1)
          for (int q = 0; q < KDmn::dmn_size(); ++q)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1)
                G4_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex) = G4_val_;
  }

  static void SetUpTestCase() {
    const std::array<double, 4> cluster_basis{1., 0., 0., 1.};
    const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 2}};
    phys::domains::cluster_domain_initializer<RDmn>::execute(cluster_basis.data(),
                                                             cluster_superbasis);
  }

  const parallel::NoConcurrency concurrency_;

  const std::complex<double> G_val_;
  const std::complex<double> G4_val_;

  const std::complex<double> G4_over_GGGG;
  const std::complex<double> one_over_GG;

  SpGreensFunction G_;
  TpGreensFunction G4_;
  TpGreensFunction Gamma_;
};

TEST_F(ComputeReducibleTpVertexTest, ParticleHoleLongitudinalUpDown) {
  phys::computeReducibleTpVertex(beta_, concurrency_, G_, G4_,
                                 phys::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN, Gamma_);

  for (int w_ex = 0; w_ex < FreqExchangeDmn::dmn_size(); ++w_ex)
    for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2)
      for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1)
        for (int q = 0; q < KDmn::dmn_size(); ++q)
          for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
            for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
              if (q == 0 && w_ex == 0)
                EXPECT_EQ(G4_over_GGGG - one_over_GG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
              else
                EXPECT_EQ(G4_over_GGGG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
            }
}

TEST_F(ComputeReducibleTpVertexTest, ParticleHoleLongitudinalUpUp) {
  phys::computeReducibleTpVertex(beta_, concurrency_, G_, G4_,
                                 phys::PARTICLE_HOLE_LONGITUDINAL_UP_UP, Gamma_);

  for (int w_ex = 0; w_ex < FreqExchangeDmn::dmn_size(); ++w_ex)
    for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2)
      for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1)
        for (int q = 0; q < KDmn::dmn_size(); ++q)
          for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
            for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
              if (q == 0 && w_ex == 0 && (k1 != k2 || w1 != w2))
                EXPECT_EQ(G4_over_GGGG - one_over_GG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
              else if ((q != 0 || w_ex != 0) && k1 == k2 && w1 == w2)
                EXPECT_EQ(G4_over_GGGG + one_over_GG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
              else {
                EXPECT_NEAR(G4_over_GGGG.real(), Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex).real(),
                            std::numeric_limits<double>::epsilon());
                EXPECT_NEAR(G4_over_GGGG.imag(), Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex).imag(),
                            std::numeric_limits<double>::epsilon());
              }
            }
}

TEST_F(ComputeReducibleTpVertexTest, ParticleHoleTransverse) {
  phys::computeReducibleTpVertex(beta_, concurrency_, G_, G4_, phys::PARTICLE_HOLE_TRANSVERSE,
                                 Gamma_);

  for (int w_ex = 0; w_ex < FreqExchangeDmn::dmn_size(); ++w_ex)
    for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2)
      for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1)
        for (int q = 0; q < KDmn::dmn_size(); ++q)
          for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
            for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
              if (k1 == k2 && w1 == w2)
                EXPECT_EQ(G4_over_GGGG + one_over_GG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
              else
                EXPECT_EQ(G4_over_GGGG, Gamma_(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex));
            }
}
