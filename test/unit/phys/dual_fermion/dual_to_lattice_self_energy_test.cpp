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
#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

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

  using Complex = DualToLatticeSelfEnergyType::Complex;

  using BandDmn = DualToLatticeSelfEnergyType::BandDmn;
  using SpinDmn = DualToLatticeSelfEnergyType::SpinDmn;

  using SpFreqDmn = DualToLatticeSelfEnergyType::SpFreqDmn;

  using KClusterDmn = DualToLatticeSelfEnergyType::KClusterDmn;
  using RClusterDmn = func::dmn_0<KClusterDmn::parameter_type::dual_type>;

  using KLatticeDmn = DualToLatticeSelfEnergyType::KLatticeDmn;
  using RLatticeDmn = func::dmn_0<KLatticeDmn::parameter_type::dual_type>;

  using KSuperlatticeDmn = DualToLatticeSelfEnergyType::KSuperlatticeDmn;
  using RSuperlatticeDmn = func::dmn_0<KSuperlatticeDmn::parameter_type::dual_type>;

  using SpClusterGF = DualToLatticeSelfEnergyType::SpClusterGF;
  using SpLatticeGF = DualToLatticeSelfEnergyType::SpLatticeGF;
  using DualGF = DualToLatticeSelfEnergyType::DualGF;

  DualToLatticeSelfEnergyTest()
      : concurrency_(0, nullptr),
        dual_to_lattice_comp_(concurrency_, G_cluster_, Sigma_cluster_, Sigma_dual_,
                              Sigma_lattice_) {}

  static void SetUpTestCase() {
    const int num_bands = 1;
    BandDmn::parameter_type::initialize(parameters_, num_bands, std::vector<int>(num_bands),
                                        std::vector<std::vector<double>>(num_bands));

    SpFreqDmn::parameter_type::initialize(parameters_);

    const std::array<double, 4> cluster_basis{1., 0., 0., 1.};
    const std::vector<std::vector<int>> cluster_superbasis{{2, 0}, {0, 1}};
    phys::domains::cluster_domain_initializer<
        func::dmn_0<KClusterDmn::parameter_type::dual_type>>::execute(cluster_basis.data(),
                                                                      cluster_superbasis);

    const std::array<double, 4> superlattice_basis{2., 0., 0., 1.};
    const int superlattice_size = 2;
    const std::vector<std::vector<int>> superlattice_superbasis{{superlattice_size, 0},
                                                                {0, superlattice_size}};
    phys::domains::cluster_domain_initializer<
        func::dmn_0<KSuperlatticeDmn::parameter_type::dual_type>>::execute(superlattice_basis.data(),
                                                                           superlattice_superbasis);

    std::vector<std::vector<int>> lattice_superbasis(cluster_superbasis);
    for (int i = 0; i < dimension_; ++i) {
      for (int d = 0; d < dimension_; ++d) {
        lattice_superbasis[i][d] *= superlattice_size;
      }
    }
    phys::domains::cluster_domain_initializer<
        func::dmn_0<KLatticeDmn::parameter_type::dual_type>>::execute(cluster_basis.data(),
                                                                      lattice_superbasis);
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

TEST_F(DualToLatticeSelfEnergyTest, ComputeNonDiagonalLatticeSelfEnergy) {
  // Prepare input.
  Complex Sigma_dual_00(1., 1.);
  Complex Sigma_dual_01(2., 4.);
  Complex Sigma_dual_10(3., 9.);
  Complex Sigma_dual_11(4., 16.);
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      Sigma_dual_(0, 0, k_tilde, w) = Sigma_dual_00;
      Sigma_dual_(0, 1, k_tilde, w) = Sigma_dual_01;
      Sigma_dual_(1, 0, k_tilde, w) = Sigma_dual_10;
      Sigma_dual_(1, 1, k_tilde, w) = Sigma_dual_11;
    }
  }

  Complex Sigma_cluster_0(3.14, 1.59);
  Complex Sigma_cluster_1(2.65, 3.59);
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int s = 0; s < SpinDmn::dmn_size(); ++s) {
      Sigma_cluster_(0, s, 0, s, 0, w) = Sigma_cluster_0;
      Sigma_cluster_(0, s, 0, s, 1, w) = Sigma_cluster_1;
    }
  }

  Complex G_cluster_0(4.2, 2.4);
  Complex G_cluster_1(3.7, 7.3);
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int s = 0; s < SpinDmn::dmn_size(); ++s) {
      G_cluster_(0, s, 0, s, 0, w) = G_cluster_0;
      G_cluster_(0, s, 0, s, 1, w) = G_cluster_1;
    }
  }

  // Compute expected result.
  const auto det_A_inv =
      1. / ((1. + Sigma_dual_00 * G_cluster_0) * (1. + Sigma_dual_11 * G_cluster_1) -
            Sigma_dual_01 * G_cluster_1 * Sigma_dual_10 * G_cluster_0);

  const auto Sigma_lattice_nondiag_K_00 =
      Sigma_cluster_0 + det_A_inv * (Sigma_dual_00 + Sigma_dual_00 * Sigma_dual_11 * G_cluster_1 -
                                     Sigma_dual_01 * Sigma_dual_10 * G_cluster_1);

  const auto Sigma_lattice_nondiag_K_01 = det_A_inv * Sigma_dual_01;

  const auto Sigma_lattice_nondiag_K_10 = det_A_inv * Sigma_dual_10;

  const auto Sigma_lattice_nondiag_K_11 =
      Sigma_cluster_1 + det_A_inv * (Sigma_dual_11 + Sigma_dual_00 * Sigma_dual_11 * G_cluster_0 -
                                     Sigma_dual_01 * Sigma_dual_10 * G_cluster_0);

  dual_to_lattice_comp_.computeNonDiagonalLatticeSelfEnergy();
  const auto& Sigma_lattice_nondiag_K = dual_to_lattice_comp_.getNonDiagonalLatticeSelfEnergy();

  // Check results.
  const double tol = std::numeric_limits<double>::epsilon();

  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      EXPECT_NEAR(Sigma_lattice_nondiag_K_00.real(),
                  Sigma_lattice_nondiag_K(0, 0, k_tilde, w).real(), tol);
      EXPECT_NEAR(Sigma_lattice_nondiag_K_00.imag(),
                  Sigma_lattice_nondiag_K(0, 0, k_tilde, w).imag(), tol);

      EXPECT_NEAR(Sigma_lattice_nondiag_K_01.real(),
                  Sigma_lattice_nondiag_K(0, 1, k_tilde, w).real(), tol);
      EXPECT_NEAR(Sigma_lattice_nondiag_K_01.imag(),
                  Sigma_lattice_nondiag_K(0, 1, k_tilde, w).imag(), tol);

      EXPECT_NEAR(Sigma_lattice_nondiag_K_10.real(),
                  Sigma_lattice_nondiag_K(1, 0, k_tilde, w).real(), tol);
      EXPECT_NEAR(Sigma_lattice_nondiag_K_10.imag(),
                  Sigma_lattice_nondiag_K(1, 0, k_tilde, w).imag(), tol);

      EXPECT_NEAR(Sigma_lattice_nondiag_K_11.real(),
                  Sigma_lattice_nondiag_K(1, 1, k_tilde, w).real(), tol);
      EXPECT_NEAR(Sigma_lattice_nondiag_K_11.imag(),
                  Sigma_lattice_nondiag_K(1, 1, k_tilde, w).imag(), tol);
    }
  }
}

TEST_F(DualToLatticeSelfEnergyTest, ComputeDiagonalLatticeSelfEnergy) {
  // Prepare diagonal (in K) input.
  DualGF Sigma_lattice_nondiag;
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      for (int K = 0; K < KClusterDmn::dmn_size(); ++K) {
        const Complex value(K + std::sin(w), k_tilde - std::cos(w));
        Sigma_lattice_nondiag(K, K, k_tilde, w) = value;
      }
    }
  }
  dual_to_lattice_comp_.setNonDiagonalLatticeSelfEnergy(Sigma_lattice_nondiag);

  // Compute expected result.
  // - Copy diagonal into single argument function.
  func::function<Complex, func::dmn_variadic<KClusterDmn, KSuperlatticeDmn, SpFreqDmn>>
      Sigma_lattice_K_k_tilde;
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int k_tilde = 0; k_tilde < KSuperlatticeDmn::dmn_size(); ++k_tilde) {
      for (int K = 0; K < KClusterDmn::dmn_size(); ++K) {
        Sigma_lattice_K_k_tilde(K, k_tilde, w) = Sigma_lattice_nondiag(K, K, k_tilde, w);
      }
    }
  }

  // - Separate momentum to real space FTs for cluster and superlattice arguments.
  func::function<Complex, func::dmn_variadic<RClusterDmn, KSuperlatticeDmn, SpFreqDmn>>
      Sigma_lattice_R_k_tilde;
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(Sigma_lattice_K_k_tilde,
                                                                        Sigma_lattice_R_k_tilde);
  func::function<Complex, func::dmn_variadic<RClusterDmn, RSuperlatticeDmn, SpFreqDmn>>
      Sigma_lattice_R_r_tilde;
  math::transform::FunctionTransform<KSuperlatticeDmn, RSuperlatticeDmn>::execute(
      Sigma_lattice_R_k_tilde, Sigma_lattice_R_r_tilde);

  // - Combine cluster and superlattice arguments into single lattice argument and add band and spin
  //   dependecy.
  func::function<Complex,
                 func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, RLatticeDmn, SpFreqDmn>>
      Sigma_lattice_r;
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w) {
    for (int r_tilde = 0; r_tilde < RSuperlatticeDmn::dmn_size(); ++r_tilde) {
      const auto& r_tilde_vec = RSuperlatticeDmn::get_elements()[r_tilde];

      for (int R = 0; R < RClusterDmn::dmn_size(); ++R) {
        const auto& R_vec = RClusterDmn::get_elements()[R];

        const auto R_plus_r_tilde = math::util::add(R_vec, r_tilde_vec);

        const auto r = phys::domains::cluster_operations::index(
            phys::domains::cluster_operations::find_closest_cluster_vector(
                R_plus_r_tilde, RLatticeDmn::get_elements(),
                RLatticeDmn::parameter_type::get_super_basis_vectors()),
            RLatticeDmn::get_elements(), RLatticeDmn::parameter_type::SHAPE);

        Sigma_lattice_r(0, 0, 0, 0, r, w) = Sigma_lattice_R_r_tilde(R, r_tilde, w);
        Sigma_lattice_r(0, 1, 0, 1, r, w) = Sigma_lattice_R_r_tilde(R, r_tilde, w);
      }
    }
  }

  // - Lattice FT to momentum space.
  SpLatticeGF Sigma_lattice_expected;
  math::transform::FunctionTransform<RLatticeDmn, KLatticeDmn>::execute(Sigma_lattice_r,
                                                                        Sigma_lattice_expected);

  dual_to_lattice_comp_.computeDiagonalLatticeSelfEnergy();

  // Check results.
  for (int w = 0; w < SpFreqDmn::dmn_size(); ++w)
    for (int k = 0; k < KLatticeDmn::dmn_size(); ++k)
      for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
        for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
          for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
            for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1) {
              EXPECT_DOUBLE_EQ(Sigma_lattice_expected(b1, s1, b2, s2, k, w).real(),
                               Sigma_lattice_(b1, s1, b2, s2, k, w).real());
              EXPECT_DOUBLE_EQ(Sigma_lattice_expected(b1, s1, b2, s2, k, w).imag(),
                               Sigma_lattice_(b1, s1, b2, s2, k, w).imag());
            }
}
