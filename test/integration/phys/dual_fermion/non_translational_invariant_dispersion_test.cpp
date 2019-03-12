// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file computes and tests the non-translational invariant dispersion required for DCA + DF.
// The test model is a 2x2 square lattice with a 4x4 superlattice.

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

TEST(NonTranslationalInvariantDispersionTest, FourSiteSquareLattice) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;
  params.set_t(1.);
  params.set_t_prime(0.);

  constexpr int dimension = 2;

  // 2x2 square cluster
  const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 2}};
  // cluster superbasis = superlattice basis
  const std::array<double, 4> superlattice_basis{2., 0., 0., 2.};

  // We build the superlattice out of 4x4=16 clusters.
  const int superlattice_size = 4;

  const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                              {0, superlattice_size}};

  std::vector<std::vector<int>> lattice_superbasis(cluster_superbasis);
  for (int i = 0; i < dimension; ++i) {
    for (int d = 0; d < dimension; ++d) {
      lattice_superbasis[i][d] *= superlattice_size;
    }
  }

  // Setup cluster.
  using RClusterDmn = phys::ClusterDomainAliases<dimension>::RClusterDmn;
  using KClusterDmn = phys::ClusterDomainAliases<dimension>::KClusterDmn;
  phys::domains::cluster_domain_initializer<RClusterDmn>::execute(Lattice::initialize_r_DCA_basis(),
                                                                  cluster_superbasis);
  // math::util::print(RClusterDmn::get_elements());
  // math::util::print(KClusterDmn::get_elements());

  // Setup lattice.
  using RLatticeDmn = phys::ClusterDomainAliases<dimension>::RSpHostDmn;
  phys::domains::cluster_domain_initializer<RLatticeDmn>::execute(Lattice::initialize_r_DCA_basis(),
                                                                  lattice_superbasis);

  // Setup superlattice.
  using RSuperlatticeDmn = phys::ClusterDomainAliases<dimension>::RSpSuperlatticeDmn;
  using KSuperlatticeDmn = phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;
  phys::domains::cluster_domain_initializer<RSuperlatticeDmn>::execute(superlattice_basis.data(),
                                                                       superlattice_superbasis);

  //////////////////////////////////////////////////////////////////////////////
  // Initialize hopping matrix in real space.
  func::function<std::complex<double>, func::dmn_variadic<RClusterDmn, RClusterDmn, RSuperlatticeDmn>>
      t_IJ_d_tilde;
  Lattice::initializeNonTranslationalInvariantHoppingMatrix<RLatticeDmn>(params, t_IJ_d_tilde);

  //////////////////////////////////////////////////////////////////////////////
  // Fourier transform from superlattice real space (d_tilde) to momentum space (k_tilde).
  func::function<std::complex<double>, func::dmn_variadic<RClusterDmn, RClusterDmn, KSuperlatticeDmn>>
      t_IJ_k_tilde;
  math::transform::FunctionTransform<RSuperlatticeDmn, KSuperlatticeDmn>::execute(t_IJ_d_tilde,
                                                                                  t_IJ_k_tilde);

  // Check t(I, J, k_tilde).
  const std::complex<double> i(0., 1.);
  const double t = params.get_t();

  for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
    const auto k_tilde_x = KSuperlatticeDmn::get_elements()[k_tilde][0];
    const auto k_tilde_y = KSuperlatticeDmn::get_elements()[k_tilde][1];

    const auto exp_2i_kx = std::exp(2. * i * k_tilde_x);
    const auto exp_min_2i_kx = std::exp(-2. * i * k_tilde_x);
    const auto exp_2i_ky = std::exp(2. * i * k_tilde_y);
    const auto exp_min_2i_ky = std::exp(-2. * i * k_tilde_y);

    EXPECT_NEAR(0., t_IJ_k_tilde(0, 0, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_ky).real(), t_IJ_k_tilde(0, 1, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_kx).real(), t_IJ_k_tilde(0, 2, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(0, 3, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(0., t_IJ_k_tilde(0, 0, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_ky).imag(), t_IJ_k_tilde(0, 1, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_kx).imag(), t_IJ_k_tilde(0, 2, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(0, 3, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(t * (1. + exp_min_2i_ky).real(), t_IJ_k_tilde(1, 0, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 1, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 2, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_kx).real(), t_IJ_k_tilde(1, 3, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(t * (1. + exp_min_2i_ky).imag(), t_IJ_k_tilde(1, 0, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 1, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 2, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_kx).imag(), t_IJ_k_tilde(1, 3, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(t * (1. + exp_min_2i_kx).real(), t_IJ_k_tilde(2, 0, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(2, 1, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(2, 2, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_ky).real(), t_IJ_k_tilde(2, 3, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(t * (1. + exp_min_2i_kx).imag(), t_IJ_k_tilde(2, 0, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(2, 1, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(2, 2, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_2i_ky).imag(), t_IJ_k_tilde(2, 3, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(0., t_IJ_k_tilde(3, 0, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_min_2i_kx).real(), t_IJ_k_tilde(3, 1, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_min_2i_ky).real(), t_IJ_k_tilde(3, 2, k_tilde).real(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(3, 3, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());

    EXPECT_NEAR(0., t_IJ_k_tilde(3, 0, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_min_2i_kx).imag(), t_IJ_k_tilde(3, 1, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(t * (1. + exp_min_2i_ky).imag(), t_IJ_k_tilde(3, 2, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(3, 3, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
  }

  //////////////////////////////////////////////////////////////////////////////
  // Fourier transform from cluster real space (I, J) to momentum space (K, K_prime).
  math::transform::SpaceTransform2D<RClusterDmn, KClusterDmn, double>::execute(t_IJ_k_tilde);

  // Check t(K, K_prime, k_tilde).
  for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
    const auto k_tilde_x = KSuperlatticeDmn::get_elements()[k_tilde][0];
    const auto k_tilde_y = KSuperlatticeDmn::get_elements()[k_tilde][1];

    const auto one_over_Nc = 1. / KClusterDmn::dmn_size();

    // t(K=0, K_prime=0, k_tilde).
    EXPECT_NEAR(one_over_Nc * 4. * t * (2 + std::cos(2. * k_tilde_x) + std::cos(2. * k_tilde_y)),
                t_IJ_k_tilde(0, 0, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(0, 0, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());

    // t(K=0, K_prime=1, k_tilde).
    EXPECT_NEAR(0., t_IJ_k_tilde(0, 1, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(-one_over_Nc * 4. * t * std::sin(2. * k_tilde_y),
                t_IJ_k_tilde(0, 1, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());

    // t(K=1, K_prime=0, k_tilde).
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 0, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(one_over_Nc * 4. * t * std::sin(2. * k_tilde_y), t_IJ_k_tilde(1, 0, k_tilde).imag(),
                100 * std::numeric_limits<double>::epsilon());

    // t(K=1, K_prime=1, k_tilde).
    EXPECT_NEAR(one_over_Nc * 4. * t * (std::cos(2. * k_tilde_x) - std::cos(2. * k_tilde_y)),
                t_IJ_k_tilde(1, 1, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(1, 1, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());

    // t(K=3, K_prime=3, k_tilde).
    EXPECT_NEAR(-one_over_Nc * 4. * t * (2 + std::cos(2. * k_tilde_x) + std::cos(2. * k_tilde_y)),
                t_IJ_k_tilde(3, 3, k_tilde).real(), 100 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0., t_IJ_k_tilde(3, 3, k_tilde).imag(), 100 * std::numeric_limits<double>::epsilon());
  }
}
