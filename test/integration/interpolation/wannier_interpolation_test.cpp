// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// The file provides tests and examples of Wannier interpolation between different momentum space
// grids.

#include <array>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

class WannierInterpolationTest : public ::testing::Test {
protected:
  // Cluster
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

  // Host
  using RHostType = dca::phys::domains::cluster_domain<double, 2, dca::phys::domains::LATTICE_SP,
                                                       dca::phys::domains::REAL_SPACE,
                                                       dca::phys::domains::BRILLOUIN_ZONE>;
  using RHostDmn = dca::func::dmn_0<RHostType>;

  using KHostType = dca::phys::domains::cluster_domain<double, 2, dca::phys::domains::LATTICE_SP,
                                                       dca::phys::domains::MOMENTUM_SPACE,
                                                       dca::phys::domains::BRILLOUIN_ZONE>;
  using KHostDmn = dca::func::dmn_0<KHostType>;

  static void SetUpTestCase() {
    std::array<double, 4> basis{{1., 0., 0., 1.}};  // Lattice basis: [1, 0], [0, 1].

    std::vector<std::vector<int>> cluster{{4, 0}, {0, 4}};  // # cluster sites = 16.
    dca::phys::domains::cluster_domain_initializer<RClusterDmn>::execute(basis.data(), cluster);

    CenteredRClusterType::initialize();

    std::vector<std::vector<int>> host{{8, 0}, {0, 8}};  // # host sites = 64.
    dca::phys::domains::cluster_domain_initializer<RHostDmn>::execute(basis.data(), host);
  }
};

TEST_F(WannierInterpolationTest, ClusterToCluster) {
  dca::func::function<std::complex<double>, KClusterDmn> f_k_cluster("K-cluster");
  for (std::size_t i = 0; i < f_k_cluster.size(); ++i)
    f_k_cluster(i) = std::complex<double>(i, i * i);
  // f_k_cluster.print_elements();

  dca::func::function<std::complex<double>, CenteredRClusterDmn> f_centered_r_cluster(
      "Centered-R-cluster");
  dca::func::function<std::complex<double>, KClusterDmn> f_k_cluster_transformed(
      "K-cluster-transformed");

  dca::math::transform::FunctionTransform<KClusterDmn, CenteredRClusterDmn>::execute(
      f_k_cluster, f_centered_r_cluster);
  // f_centered_r_cluster.print_elements();

  dca::math::transform::FunctionTransform<CenteredRClusterDmn, KClusterDmn>::execute(
      f_centered_r_cluster, f_k_cluster_transformed);
  // f_k_cluster_transformed.print_elements();

  for (std::size_t i = 0; i < f_k_cluster.size(); ++i) {
    EXPECT_NEAR(f_k_cluster(i).real(), f_k_cluster_transformed(i).real(),
                3000 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(f_k_cluster(i).imag(), f_k_cluster_transformed(i).imag(),
                3000 * std::numeric_limits<double>::epsilon());
  }
}

TEST_F(WannierInterpolationTest, ClusterToHost) {
  dca::func::function<std::complex<double>, KClusterDmn> f_k_cluster("K-cluster");
  for (std::size_t i = 0; i < f_k_cluster.size(); ++i)
    f_k_cluster(i) = std::complex<double>(i, i * i);
  // f_k_cluster.print_elements();

  dca::func::function<std::complex<double>, CenteredRClusterDmn> f_centered_r_cluster(
      "Centered-R-cluster");
  dca::func::function<std::complex<double>, KHostDmn> f_k_host("K-host");

  dca::math::transform::FunctionTransform<KClusterDmn, CenteredRClusterDmn>::execute(
      f_k_cluster, f_centered_r_cluster);
  // f_centered_r_cluster.print_elements();

  dca::math::transform::FunctionTransform<CenteredRClusterDmn, KHostDmn>::execute(
      f_centered_r_cluster, f_k_host);
  // f_k_host.print_elements();

  for (std::size_t cluster_ind = 0; cluster_ind < KClusterDmn::dmn_size(); ++cluster_ind) {
    const auto& k_cluster = KClusterType::get_elements()[cluster_ind];

    const std::size_t host_ind = dca::phys::domains::cluster_operations::index(
        k_cluster, KHostType::get_elements(), KHostType::SHAPE);

    EXPECT_NEAR(f_k_cluster(cluster_ind).real(), f_k_host(host_ind).real(),
                3000 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(f_k_cluster(cluster_ind).imag(), f_k_host(host_ind).imag(),
                3000 * std::numeric_limits<double>::epsilon());
  }
}
