// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dual_to_lattice_self_energy.hpp.

#include "dca/phys/dual_fermion/dual_to_lattice_self_energy.hpp"

#include <array>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"

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
    return 4;
  }
};
}  // namespace testing

class DualToLatticeSelfEnergyTest : public ::testing::Test {
protected:
  static constexpr int dimension_ = 2;

  using DualToLatticeSelfEnergyType =
      phys::df::DualToLatticeSelfEnergy<double, parallel::NoConcurrency, dimension_>;

  using BandDmn = DualToLatticeSelfEnergyType::BandDmn;
  using SpinDmn = DualToLatticeSelfEnergyType::SpinDmn;

  using SpFreqDmn = DualToLatticeSelfEnergyType::SpFreqDmn;

  using KClusterDmn = DualToLatticeSelfEnergyType::KClusterDmn;
  using KLatticeDmn = DualToLatticeSelfEnergyType::KLatticeDmn;
  using KSuperlatticeDmn = DualToLatticeSelfEnergyType::KSuperlatticeDmn;

  using SpClusterGF = DualToLatticeSelfEnergyType::SpClusterGF;
  using SpLatticeGF = DualToLatticeSelfEnergyType::SpLatticeGF;
  using DualGF = DualToLatticeSelfEnergyType::DualGF;

  DualToLatticeSelfEnergyTest()
      : concurrency_(0, nullptr),
        dual_to_lattice_comp_(concurrency_, G_cluster_, Sigma_cluster_, Sigma_dual_, Sigma_lattice_) {
  }

  static void SetUpTestCase() {
    const int num_bands = 1;
    BandDmn::parameter_type::initialize(parameters_, num_bands, std::vector<int>(num_bands),
                                        std::vector<std::vector<double>>(num_bands));

    SpFreqDmn::parameter_type::initialize(parameters_);

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

    std::vector<std::vector<int>> lattice_superbasis(cluster_superbasis);
    for (int i = 0; i < dimension_; ++i) {
      for (int d = 0; d < dimension_; ++d) {
        lattice_superbasis[i][d] *= superlattice_size;
      }
    }
    phys::domains::cluster_domain_initializer<func::dmn_0<KLatticeDmn::parameter_type::dual_type>>::execute(
        cluster_basis.data(), lattice_superbasis);
  }

  static const testing::MockParameters parameters_;

  const parallel::NoConcurrency concurrency_;

  SpClusterGF G_cluster_;
  SpClusterGF Sigma_cluster_;
  DualGF Sigma_dual_;
  SpLatticeGF Sigma_lattice_;

  DualToLatticeSelfEnergyType dual_to_lattice_comp_;
};

const testing::MockParameters DualToLatticeSelfEnergyTest::parameters_;

TEST_F(DualToLatticeSelfEnergyTest, ComputeNonInvariantLatticeSelfEnergy) {}

TEST_F(DualToLatticeSelfEnergyTest, FourierTransformToRealSpace) {}

TEST_F(DualToLatticeSelfEnergyTest, FourierTransformToMomentumSpace) {}
