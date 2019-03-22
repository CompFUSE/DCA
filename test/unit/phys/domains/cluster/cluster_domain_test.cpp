// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file illustrates relation between and use of various cluster domains.

#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

#include <array>
#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

using namespace dca;

TEST(ClusterDomainTest, ClusterLatticeSuperlattice) {
  constexpr int dimension = 2;

  // Basis spanning the lattice.
  const std::array<double, 4> basis{1., 0., 0., 1.};

  // Superbasis defining the cluster (cluster superbasis = superlattice basis).
  const std::vector<std::vector<int>> cluster_superbasis{{2, 2}, {2, -2}};
  // Need a std::array to pass to cluster_domain_initializer::execute.
  const std::array<double, 4> superlattice_basis{2., 2., 2., -2.};

  // Linear size of the superlattice.
  const int superlattice_size = 4;
  // Expressed as superbasis vectors:
  const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                              {0, superlattice_size}};

  // The superbasis defining the lattice is obtained by scaling the cluster superbasis.
  std::vector<std::vector<int>> lattice_superbasis(cluster_superbasis);
  for (int i = 0; i < dimension; ++i) {
    for (int d = 0; d < dimension; ++d) {
      lattice_superbasis[i][d] *= superlattice_size;
    }
  }

  // Setup cluster.
  using RClusterDmn = phys::ClusterDomainAliases<dimension>::RClusterDmn;
  using KClusterDmn = phys::ClusterDomainAliases<dimension>::KClusterDmn;

  phys::domains::cluster_domain_initializer<RClusterDmn>::execute(basis.data(), cluster_superbasis);

  const std::vector<std::vector<double>>& R = RClusterDmn::get_elements();
  const std::vector<std::vector<double>>& K = KClusterDmn::get_elements();
  // math::util::print(R);
  // math::util::print(K);

  // Setup lattice.
  using RLatticeDmn = phys::ClusterDomainAliases<dimension>::RSpHostDmn;
  using KLatticeDmn = phys::ClusterDomainAliases<dimension>::KSpHostDmn;

  phys::domains::cluster_domain_initializer<RLatticeDmn>::execute(basis.data(), lattice_superbasis);

  const std::vector<std::vector<double>>& r = RLatticeDmn::get_elements();
  const std::vector<std::vector<double>>& k = KLatticeDmn::get_elements();
  // math::util::print(r);
  // math::util::print(k);

  // Setup superlattice.
  using RSuperlatticeDmn = phys::ClusterDomainAliases<dimension>::RSpSuperlatticeDmn;
  using KSuperlatticeDmn = phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;

  phys::domains::cluster_domain_initializer<RSuperlatticeDmn>::execute(superlattice_basis.data(),
                                                                       superlattice_superbasis);

  const std::vector<std::vector<double>>& r_tilde = RSuperlatticeDmn::get_elements();
  const std::vector<std::vector<double>>& k_tilde = KSuperlatticeDmn::get_elements();
  // math::util::print(r_tilde);
  // math::util::print(k_tilde);

  // Setup centered reciprocal superlattice.
  using CenteredKSuperlatticeType =
      phys::ClusterDomainAliases<dimension>::CenteredKSpSuperlatticeType;

  CenteredKSuperlatticeType::initialize();

  const std::vector<std::vector<double>>& k_tilde_centered =
      CenteredKSuperlatticeType::get_elements();
  // const std::vector<double>& weights = CenteredKSuperlatticeType::get_weights();
  // math::util::print(k_tilde_centered);
  // math::util::print(weights);

  // Check the domain sizes.
  EXPECT_EQ(KLatticeDmn::dmn_size(), KClusterDmn::dmn_size() * KSuperlatticeDmn::dmn_size());

  // Check that real space cluster vectors are also real space lattice vectors.
  const double tolerance = 1.e-6;
  for (const auto& R_vec : R) {
    phys::domains::cluster_operations::find_closest_cluster_vector(
        R_vec, r, RLatticeDmn::parameter_type::get_super_basis_vectors(), tolerance);
  }

  // Check that reciprocal cluster vectors are also reciprocal lattice vectors.
  for (const auto& K_vec : K) {
    phys::domains::cluster_operations::find_closest_cluster_vector(
        K_vec, k, KLatticeDmn::parameter_type::get_super_basis_vectors(), tolerance);
  }

  // Check that real space superlattice vectors are also real space lattice vectors.
  for (const auto& r_tilde_vec : r_tilde) {
    phys::domains::cluster_operations::find_closest_cluster_vector(
        r_tilde_vec, r, RLatticeDmn::parameter_type::get_super_basis_vectors(), tolerance);
  }

  // Check that reciprocal superlattice vectors are also reciprocal lattice vectors.
  for (const auto& k_tilde_vec : k_tilde) {
    phys::domains::cluster_operations::find_closest_cluster_vector(
        k_tilde_vec, k, KLatticeDmn::parameter_type::get_super_basis_vectors(), tolerance);
  }

  // Check that centered reciprocal superlattice vectors are also reciprocal lattice vectors.
  for (const auto& k_tilde_vec : k_tilde_centered) {
    phys::domains::cluster_operations::find_closest_cluster_vector(
        k_tilde_vec, k, KLatticeDmn::parameter_type::get_super_basis_vectors(), tolerance);
  }
}
