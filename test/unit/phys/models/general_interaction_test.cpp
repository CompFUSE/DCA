// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests general_interaction.hpp.

#include "dca/phys/models/general_interaction.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

// Provide a parameters class that extends ModelParameters by the static data member
// 'lattice_dimension' and the member function 'get_interacting_orbitals', which are required by
// general_interaction.
template <typename Lattice>
class FakeParameters
    : public phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> {
public:
  constexpr static int lattice_dimension = Lattice::DIMENSION;

  void set_interacting_orbitals(const std::vector<int>& interacting_orbitals) {
    interacting_orbitals_ = interacting_orbitals;
  }
  const std::vector<int>& get_interacting_orbitals() const {
    return interacting_orbitals_;
  }

private:
  std::vector<int> interacting_orbitals_;
};

TEST(GeneralInteractionTest, MakeCorrelatedOrbitals) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::bilayer_lattice<PointGroup>;

  using BandDmn = func::dmn_0<func::dmn<2, int>>;
  using SpinDmn = func::dmn_0<func::dmn<2, int>>;
  using BandSpinDmn = func::dmn_variadic<BandDmn, SpinDmn>;

  using CDA = phys::ClusterDomainAliases<Lattice::DIMENSION>;
  using RClusterDmn= typename CDA::RClusterDmn;

  const std::vector<std::vector<int>> DCA_cluster{{1, 0}, {0, 1}};
  phys::domains::cluster_domain_initializer<RClusterDmn>::execute(Lattice::initialize_r_DCA_basis(),
                                                            DCA_cluster);

  using ParametersType = FakeParameters<Lattice>;
  ParametersType params;
  params.set_U(4);
  params.set_V(2);
  params.set_V_prime(1);

  func::function<double, func::dmn_variadic<BandSpinDmn, BandSpinDmn, RClusterDmn>> H_interaction;
  Lattice::initialize_H_interaction(H_interaction, params);

  // Both orbitals are interacting.
  params.set_interacting_orbitals({0, 1});
  std::vector<int> correlated_orbitals =
      phys::models::general_interaction<ParametersType>::make_correlated_orbitals(params,
                                                                                  H_interaction);

  std::vector<int> correlated_orbitals_check{1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14};
  EXPECT_EQ(correlated_orbitals_check, correlated_orbitals);

  // Only second orbital (index=1) is interacting.
  params.set_interacting_orbitals({1});
  correlated_orbitals = phys::models::general_interaction<ParametersType>::make_correlated_orbitals(
      params, H_interaction);

  correlated_orbitals_check = {7, 13};
  EXPECT_EQ(correlated_orbitals_check, correlated_orbitals);
}
