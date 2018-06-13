// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides tests and examples of Fourier transformations between the real and momentum
// space clusters.

#include <array>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"

class ClusterFourierTransformTest : public ::testing::Test {
protected:
  using RClusterType = dca::phys::domains::cluster_domain<double, 2, dca::phys::domains::CLUSTER,
                                                          dca::phys::domains::REAL_SPACE,
                                                          dca::phys::domains::BRILLOUIN_ZONE>;
  using RClusterDmn = dca::func::dmn_0<RClusterType>;

  using CenteredRClusterType = dca::phys::domains::centered_cluster_domain<RClusterType>;
  using CenteredRClusterDmn = dca::func::dmn_0<CenteredRClusterType>;

  using KClusterType = dca::phys::domains::cluster_domain<double, 2, dca::phys::domains::CLUSTER,
                                                          dca::phys::domains::MOMENTUM_SPACE,
                                                          dca::phys::domains::BRILLOUIN_ZONE>;
  using KClusterDmn = dca::func::dmn_0<KClusterType>;

  static void SetUpTestCase() {
    std::array<double, 4> basis{{1., 0., 0., 1.}};           // Lattice basis: [1, 0], [0, 1].
    std::vector<std::vector<int>> cluster{{2, 4}, {4, -2}};  // Nc = 20.

    dca::phys::domains::cluster_domain_initializer<RClusterDmn>::execute(basis.data(), cluster);

    CenteredRClusterType::initialize();
  }
};

// Basic test for a Fourier transformation between the momentum space cluster and its dual space,
// the real space cluster.
TEST_F(ClusterFourierTransformTest, DualSpaces) {
  // Test function:
  //     f(\vec{K}) = 1, if \vec{K} = 0,
  //                = 0, otherwise.
  // The Fourier transform of this function is a constant,
  //    \hat{f}(\vec{R}) = 1/Nc,
  // where Nc is the number of vectors \vec{K}.
  dca::func::function<double, KClusterDmn> f_k_cluster("K-cluster");
  f_k_cluster(0) = 1.;
  // f_k_cluster.print_elements();

  dca::func::function<double, RClusterDmn> f_r_cluster("R-cluster");
  dca::math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(f_k_cluster,
                                                                             f_r_cluster);
  // f_r_cluster.print_elements();

  for (std::size_t i = 0; i < f_r_cluster.size(); ++i)
    EXPECT_NEAR(1. / KClusterDmn::dmn_size(), f_r_cluster(i),
                10 * std::numeric_limits<double>::epsilon());

  // Back transformation
  dca::func::function<double, KClusterDmn> f_k_cluster_back("K-cluster-back-transformed");
  dca::math::transform::FunctionTransform<RClusterDmn, KClusterDmn>::execute(f_r_cluster,
                                                                             f_k_cluster_back);
  // f_k_cluster_back.print_elements();

  for (std::size_t i = 0; i < f_k_cluster.size(); ++i)
    EXPECT_NEAR(f_k_cluster(i), f_k_cluster_back(i), 10 * std::numeric_limits<double>::epsilon());
}

// Tests the Fourier transformation between a momentum space cluster and the corresponding real
// space cluster centered around the origin. This Fourier transformation is relevant for the Wannier
// intepolation.
TEST_F(ClusterFourierTransformTest, CenteredRCusterDomain) {
  // Same test function as in DualSpacs test (see above).
  dca::func::function<double, KClusterDmn> f_k_cluster("K-cluster");
  f_k_cluster(0) = 1.;
  // f_k_cluster.print_elements();

  dca::func::function<double, CenteredRClusterDmn> f_centered_r_cluster("Centered-R-cluster");
  dca::math::transform::FunctionTransform<KClusterDmn, CenteredRClusterDmn>::execute(
      f_k_cluster, f_centered_r_cluster);
  // f_centered_r_cluster.print_elements();

  // The real space Fourier components are scaled with the weight of the centered real space cluster
  // vectors.
  for (std::size_t i = 0; i < f_centered_r_cluster.size(); ++i)
    EXPECT_NEAR(CenteredRClusterType::get_weights()[i] / KClusterDmn::dmn_size(),
                f_centered_r_cluster(i), 10 * std::numeric_limits<double>::epsilon());

  // Back transformation
  dca::func::function<double, KClusterDmn> f_k_cluster_back("K-cluster-back-transformed");
  dca::math::transform::FunctionTransform<CenteredRClusterDmn, KClusterDmn>::execute(
      f_centered_r_cluster, f_k_cluster_back);
  // f_k_cluster_back.print_elements();

  for (std::size_t i = 0; i < f_k_cluster.size(); ++i)
    EXPECT_NEAR(f_k_cluster(i), f_k_cluster_back(i), 10 * std::numeric_limits<double>::epsilon());
}
