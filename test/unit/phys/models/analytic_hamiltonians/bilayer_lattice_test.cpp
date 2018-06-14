// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests bilayer_lattice.hpp.

#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

TEST(BilayerLatticeTest, Initialize_H_0) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::bilayer_lattice<PointGroup>;

  using BandDmn = func::dmn<2, int>;
  using SpinDmn = func::dmn<2, int>;
  using BandSpinDmn = func::dmn_variadic<func::dmn_0<BandDmn>, func::dmn_0<SpinDmn>>;

  using KDmn = func::dmn<4, std::vector<double>>;
  KDmn::set_elements({{0., 0.}, {0., M_PI}, {M_PI, 0.}, {M_PI, M_PI}});

  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, func::dmn_0<KDmn>>> H_0;

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;
  params.set_t(1.);
  params.set_t_prime(0.5);
  params.set_t_perp(0.25);

  Lattice::initialize_H_0(params, H_0);

  // All imaginary parts should be zero.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
          for (int k = 0; k < KDmn::dmn_size(); ++k)
            EXPECT_DOUBLE_EQ(0., H_0(b1, s1, b2, s2, k).imag());

  // All matrix elements with different spin indices should be zero.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
      for (int k = 0; k < KDmn::dmn_size(); ++k) {
        EXPECT_DOUBLE_EQ(0., H_0(b1, 0, b2, 1, k).real());
        EXPECT_DOUBLE_EQ(0., H_0(b1, 1, b2, 0, k).real());
      }

  // Check nonvanishing Hamiltonian matrix elements.
  EXPECT_DOUBLE_EQ(-6., H_0(0, 0, 0, 0, 0).real());
  EXPECT_DOUBLE_EQ(-6., H_0(0, 1, 0, 1, 0).real());
  EXPECT_DOUBLE_EQ(-6., H_0(1, 0, 1, 0, 0).real());
  EXPECT_DOUBLE_EQ(-6., H_0(1, 1, 1, 1, 0).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 0, 1, 0, 0).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 1, 1, 1, 0).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 0, 0, 0, 0).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 1, 0, 1, 0).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 0, 1, 0, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 1, 1, 1, 1).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 0, 1, 0, 1).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 1, 1, 1, 1).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 0, 0, 0, 1).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 1, 0, 1, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 0, 1, 0, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 1, 1, 1, 2).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 0, 1, 0, 2).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 1, 1, 1, 2).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 0, 0, 0, 2).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 1, 0, 1, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 3).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 3).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 0, 1, 0, 3).real());
  EXPECT_DOUBLE_EQ(2., H_0(1, 1, 1, 1, 3).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 0, 1, 0, 3).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(0, 1, 1, 1, 3).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 0, 0, 0, 3).real());
  EXPECT_DOUBLE_EQ(-0.25, H_0(1, 1, 0, 1, 3).real());
}

TEST(BilayerLatticeTest, Initialize_H_interaction) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::bilayer_lattice<PointGroup>;

  using BandDmn = func::dmn_0<func::dmn<2, int>>;
  using SpinDmn = func::dmn_0<func::dmn<2, int>>;
  using BandSpinDmn = func::dmn_variadic<BandDmn, SpinDmn>;

  using CDA = phys::ClusterDomainAliases<Lattice::DIMENSION>;
  using RClusterType = typename CDA::RClusterType;
  using RClusterDmn = typename CDA::RClusterDmn;

  const std::vector<std::vector<int>> DCA_cluster{{-2, 0}, {0, 2}};
  phys::domains::cluster_domain_initializer<RClusterDmn>::execute(Lattice::initialize_r_DCA_basis(),
                                                            DCA_cluster);

  // Get index of origin and check it.
  const int origin = RClusterType::origin_index();
  ASSERT_DOUBLE_EQ(0., RClusterDmn::get_elements()[origin][0]);
  ASSERT_DOUBLE_EQ(0., RClusterDmn::get_elements()[origin][1]);

  func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RClusterDmn>> H_interaction;
  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;

  // Check same band, opposite spin interaction.
  params.set_U(4);
  params.set_V(0);
  params.set_V_prime(0);

  Lattice::initialize_H_interaction(H_interaction, params);

  for (int r = 0; r < RClusterDmn::dmn_size(); ++r)
    for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
          for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
            if (r == origin && b1 == b2 && s1 != s2)
              EXPECT_DOUBLE_EQ(4., H_interaction(b1, s1, b2, s2, r));
            else
              EXPECT_DOUBLE_EQ(0., H_interaction(b1, s1, b2, s2, r));

  // Check different band, opposite spin interaction.
  params.set_U(0);
  params.set_V(2);
  params.set_V_prime(0);

  Lattice::initialize_H_interaction(H_interaction, params);

  for (int r = 0; r < RClusterDmn::dmn_size(); ++r)
    for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
          for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
            if (r == origin && b1 != b2 && s1 != s2)
              EXPECT_DOUBLE_EQ(2., H_interaction(b1, s1, b2, s2, r));
            else
              EXPECT_DOUBLE_EQ(0., H_interaction(b1, s1, b2, s2, r));

  // Check different band, same spin interaction.
  params.set_U(0);
  params.set_V(0);
  params.set_V_prime(1);

  Lattice::initialize_H_interaction(H_interaction, params);

  for (int r = 0; r < RClusterDmn::dmn_size(); ++r)
    for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
          for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
            if (r == origin && b1 != b2 && s1 == s2)
              EXPECT_DOUBLE_EQ(1., H_interaction(b1, s1, b2, s2, r));
            else
              EXPECT_DOUBLE_EQ(0., H_interaction(b1, s1, b2, s2, r));
}
