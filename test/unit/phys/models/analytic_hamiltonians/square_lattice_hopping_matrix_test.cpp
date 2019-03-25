// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the hopping matrix initialization routines of square_lattice.hpp.

#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"

#include <array>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

class SquareLatticeHoppingMatrixTest : public ::testing::Test {
protected:
  static constexpr int dimension_ = 2;

  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;

  using RClusterDmn = phys::ClusterDomainAliases<dimension_>::RClusterDmn;
  using RLatticeDmn = phys::ClusterDomainAliases<dimension_>::RSpHostDmn;
  using RSuperlatticeDmn = phys::ClusterDomainAliases<dimension_>::RSpSuperlatticeDmn;

  SquareLatticeHoppingMatrixTest()
      : params_(),
        t_IJ_d_tilde_(),
        superlattice_(RSuperlatticeDmn::get_elements()),
        superlattice_basis_vecs_(RSuperlatticeDmn::parameter_type::get_basis_vectors())

  {
    params_.set_t(1.);
    params_.set_t_prime(0.);
  }

  void initializeDomains(const std::vector<std::vector<int>>& cluster_superbasis,
                         const std::array<double, 4>& superlattice_basis,
                         const int superlattice_size) {
    const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                                {0, superlattice_size}};

    std::vector<std::vector<int>> lattice_superbasis(cluster_superbasis);
    for (int i = 0; i < SquareLatticeHoppingMatrixTest::dimension_; ++i) {
      for (int d = 0; d < SquareLatticeHoppingMatrixTest::dimension_; ++d) {
        lattice_superbasis[i][d] *= superlattice_size;
      }
    }

    // Setup cluster.
    phys::domains::cluster_domain_initializer<RClusterDmn>::execute(
        Lattice::initialize_r_DCA_basis(), cluster_superbasis);
    // math::util::print(RClusterDmn::get_elements());

    // Setup lattice.
    phys::domains::cluster_domain_initializer<RLatticeDmn>::execute(
        Lattice::initialize_r_DCA_basis(), lattice_superbasis);

    // Setup superlattice.
    phys::domains::cluster_domain_initializer<RSuperlatticeDmn>::execute(superlattice_basis.data(),
                                                                         superlattice_superbasis);
    // math::util::print(RSuperlatticeDmn::get_elements());

    t_IJ_d_tilde_.reset();
  }

  int countZeroElements() const {
    int num_zero = 0;

    for (int i = 0; i < t_IJ_d_tilde_.size(); ++i) {
      if (t_IJ_d_tilde_(i) == 0.)
        ++num_zero;
    }

    return num_zero;
  }

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params_;

  func::function<double, func::dmn_variadic<RClusterDmn, RClusterDmn, RSuperlatticeDmn>> t_IJ_d_tilde_;

  const std::vector<std::vector<double>>& superlattice_;
  const std::vector<std::vector<double>>& superlattice_basis_vecs_;
};

// Test on a 2x2 square cluster with a 4x4 superlattice
//
// Cluster vectors:
//  0: [0, 0]
//  1: [0, 1]
//  2: [1, 0]
//  3: [1, 1]
//
// Superlattice vectors:
//  0: [0, 0]
//  1: [0, 2]
//  2: [0, 4]
//  3: [0, 6]
//  4: [2, 0]
//  5: [2, 2]
//  6: [2, 4]
//  7: [2, 6]
//  8: [4, 0]
//  9: [4, 2]
// 10: [4, 4]
// 11: [4, 6]
// 12: [6, 0]
// 13: [6, 2]
// 14: [6, 4]
// 15: [6, 6]
TEST_F(SquareLatticeHoppingMatrixTest, NonTranslationalInvariant4siteSquare) {
  // 2x2 square cluster
  const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 2}};
  // cluster superbasis = superlattice basis
  const std::array<double, 4> superlattice_basis{2., 0., 0., 2.};

  // We build the superlattice out of 4x4=16 clusters.
  const int superlattice_size = 4;

  initializeDomains(cluster_superbasis, superlattice_basis, superlattice_size);

  Lattice::initializeNonTranslationalInvariantHoppingMatrix<RLatticeDmn>(params_, t_IJ_d_tilde_);

  // Check non-zero elements.
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 2, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 2, 4));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 1, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 1, 1));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 3, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 3, 4));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 0, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 0, 3));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 0, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 0, 12));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 3, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 3, 1));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 1, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 1, 12));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 2, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 2, 3));

  // Check the number of zero elements.
  const int num_zero = countZeroElements();
  const int num_non_zero = RClusterDmn::dmn_size() * 4;
  EXPECT_EQ(t_IJ_d_tilde_.size() - num_non_zero, num_zero);
}

// Test on a tilted 8-site cluster with 2x2 superlattice.
//
// Cluster vectors:
//  0: [0, 0]
//  1: [1,-1]
//  2: [1, 0]
//  3: [1, 1]
//  4: [2,-1]
//  5: [2, 0]
//  6: [2, 1]
//  7: [3, 0]
//
// Superlattice vectors:
//  0: [0, 0]
//  1: [2,-2]
//  2: [2, 2]
//  3: [4, 0]
TEST_F(SquareLatticeHoppingMatrixTest, NonTranslationalInvariant8siteTilted) {
  // 8-site tilted cluster
  const std::vector<std::vector<int>> cluster_superbasis{{2, 2}, {2, -2}};
  // cluster superbasis = superlattice basis
  const std::array<double, 4> superlattice_basis{2., 2., 2., -2.};

  // We build the superlattice out of 2x2=4 clusters.
  const int superlattice_size = 2;

  initializeDomains(cluster_superbasis, superlattice_basis, superlattice_size);

  Lattice::initializeNonTranslationalInvariantHoppingMatrix<RLatticeDmn>(params_, t_IJ_d_tilde_);

  // Check non-zero elements.
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 4, 1));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 2, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 6, 2));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(0, 7, 3));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 2, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 4, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 7, 2));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(1, 6, 2));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 3, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 5, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 1, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(2, 0, 0));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 7, 1));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 6, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 2, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(3, 4, 1));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(4, 5, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(4, 3, 1));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(4, 0, 1));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(4, 1, 0));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(5, 6, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(5, 7, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(5, 4, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(5, 2, 0));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(6, 0, 2));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(6, 1, 2));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(6, 5, 0));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(6, 3, 0));

  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(7, 1, 2));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(7, 0, 3));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(7, 3, 1));
  EXPECT_DOUBLE_EQ(params_.get_t(), t_IJ_d_tilde_(7, 5, 0));

  // Check the number of zero elements.
  const int num_zero = countZeroElements();
  const int num_non_zero = RClusterDmn::dmn_size() * 4;
  EXPECT_EQ(t_IJ_d_tilde_.size() - num_non_zero, num_zero);
}
