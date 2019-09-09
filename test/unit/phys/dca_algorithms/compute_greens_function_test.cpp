// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_greens_function.hpp.

#include "dca/phys/dca_algorithms/compute_greens_function.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/dca_algorithms/compute_free_greens_function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

class ComputeGreensFunctionTest : public ::testing::Test {
protected:
  using PointGroup = phys::domains::D4;
  using SpinDmn = func::dmn<2, int>;
  using MatsubaraFreqDmn = func::dmn<4, double>;  // Matsubara frequency domain with 4 elements.

  ComputeGreensFunctionTest() : concurrency_(0, nullptr) {}

  static void SetUpTestCase() {
    std::vector<double> freqs(4);
    freqs[0] = -3. * M_PI / beta;
    freqs[1] = -M_PI / beta;
    freqs[2] = M_PI / beta;
    freqs[3] = 3. * M_PI / beta;
    MatsubaraFreqDmn::set_elements(freqs);
  }

  static constexpr double beta = 1.;  // inverse temperature

  const parallel::NoConcurrency concurrency_;
};

// Test for a 2x2 square lattice (diagonal H_0) with vanishing self-energy.
TEST_F(ComputeGreensFunctionTest, SquareLatticeNoninteracting) {
  using Lattice = phys::models::square_lattice<PointGroup>;
  using OrbitalDmn = func::dmn<1, int>;  // 1 orbital
  using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<OrbitalDmn>, func::dmn_0<SpinDmn>>;

  // Momentum space domain of the 2x2 square lattice
  using KDmn = func::dmn<4, std::vector<double>>;
  const std::vector<std::vector<double>> k_vecs{{0., 0.}, {0., M_PI}, {M_PI, 0.}, {M_PI, M_PI}};
  KDmn::set_elements(k_vecs);

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;
  params.set_t(1.);

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>>>
      H0_k;
  Lattice::initialize_H_0(params, H0_k);

  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      S_k_w;  // self-energy = 0.
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G0_k_w;
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G_k_w;

  const double mu = 0.9;  // chemical potential

  // Vanishing self-energy: G = G_0.
  phys::compute_G0_k_w(H0_k, mu, 1, G0_k_w);
  phys::compute_G_k_w(H0_k, S_k_w, mu, 2, G_k_w);

  for (int i = 0; i < G_k_w.size(); ++i) {
    EXPECT_DOUBLE_EQ(G0_k_w(i).real(), G_k_w(i).real());
    EXPECT_DOUBLE_EQ(G0_k_w(i).imag(), G_k_w(i).imag());
  }
}

// Test for a single-site bilayer lattice with off-diagonal elements (in orbital-spin space) in the
// non-interacting Hamiltonian H_0.
TEST_F(ComputeGreensFunctionTest, BilayerLattice) {
  using Lattice = phys::models::bilayer_lattice<PointGroup>;
  using OrbitalDmn = func::dmn<2, int>;  // 2 orbitals
  using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<OrbitalDmn>, func::dmn_0<SpinDmn>>;

  // Momentum space domain of the single-site bilayer lattice
  using KDmn = func::dmn<1, std::vector<double>>;
  const std::vector<std::vector<double>> k_vecs{{0., 0.}};
  KDmn::set_elements(k_vecs);

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;
  params.set_t_perp(1.);
  params.set_U(4.);

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>>>
      H0_k;
  Lattice::initialize_H_0(params, H0_k);

  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      S_k_w;  // self-energy = 0.
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G0_k_w;
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G_k_w;

  const double mu = 0.9;  // chemical potential

  //
  // Vanishing self-energy: G = G_0
  //
  phys::compute_G0_k_w(H0_k, mu, 2, G0_k_w);
  phys::compute_G_k_w(H0_k, S_k_w, mu, 4, G_k_w);

  for (int i = 0; i < G_k_w.size(); ++i) {
    EXPECT_DOUBLE_EQ(G0_k_w(i).real(), G_k_w(i).real());
    EXPECT_DOUBLE_EQ(G0_k_w(i).imag(), G_k_w(i).imag());
  }

  //
  // Non-vanishing self-energy: check against ED results
  //
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G_k_w_check;

  // Initialze the non-zero elements of the self-energy and the Green's function (for later check)
  // with values from the ED solver.
  S_k_w(0, 0, 0, 0, 0, 0) = S_k_w(0, 1, 0, 1, 0, 0) = S_k_w(1, 0, 1, 0, 0, 0) =
      S_k_w(1, 1, 1, 1, 0, 0) = std::complex<double>(0.29865343066397, 0.39457553115474);
  S_k_w(0, 0, 0, 0, 0, 1) = S_k_w(0, 1, 0, 1, 0, 1) = S_k_w(1, 0, 1, 0, 0, 1) =
      S_k_w(1, 1, 1, 1, 0, 1) = std::complex<double>(0.5470450571011, 0.90951314828496);
  S_k_w(1, 0, 0, 0, 0, 0) = S_k_w(1, 1, 0, 1, 0, 0) = S_k_w(0, 0, 1, 0, 0, 0) =
      S_k_w(0, 1, 1, 1, 0, 0) = std::complex<double>(-0.01726085003162, 0.00506575883209);
  S_k_w(1, 0, 0, 0, 0, 1) = S_k_w(0, 0, 1, 0, 0, 1) = S_k_w(1, 1, 0, 1, 0, 1) =
      S_k_w(0, 1, 1, 1, 0, 1) = std::complex<double>(-0.03981884015876, 0.02265748074481);

  G_k_w_check(0, 0, 0, 0, 0, 0) = G_k_w_check(0, 1, 0, 1, 0, 0) = G_k_w_check(1, 0, 1, 0, 0, 0) =
      G_k_w_check(1, 1, 1, 1, 0, 0) = std::complex<double>(0.00600816727634, 0.10040345924141);
  G_k_w_check(0, 0, 0, 0, 0, 1) = G_k_w_check(0, 1, 0, 1, 0, 1) = G_k_w_check(1, 0, 1, 0, 0, 1) =
      G_k_w_check(1, 1, 1, 1, 0, 1) = std::complex<double>(0.01700614318379, 0.23050257305284);
  G_k_w_check(1, 0, 0, 0, 0, 0) = G_k_w_check(1, 1, 0, 1, 0, 0) = G_k_w_check(0, 0, 1, 0, 0, 0) =
      G_k_w_check(0, 1, 1, 1, 0, 0) = std::complex<double>(0.01031846228466, -0.00130614148552);
  G_k_w_check(1, 0, 0, 0, 0, 1) = G_k_w_check(1, 1, 0, 1, 0, 1) = G_k_w_check(0, 0, 1, 0, 0, 1) =
      G_k_w_check(0, 1, 1, 1, 0, 1) = std::complex<double>(0.05813527976005, -0.01071930780469);

  // Use the symmetry in Matsubara frequency to set the rest of the elements.
  for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
    for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
      for (int w = 0; w < MatsubaraFreqDmn::dmn_size() / 2; ++w) {
        S_k_w(m, n, 0, MatsubaraFreqDmn::dmn_size() - 1 - w) = std::conj(S_k_w(m, n, 0, w));
        G_k_w_check(m, n, 0, MatsubaraFreqDmn::dmn_size() - 1 - w) =
            std::conj(G_k_w_check(m, n, 0, w));
      }
    }
  }

  phys::compute_G_k_w(H0_k, S_k_w, mu, 1, G_k_w);

  const double tol = 1.e-14;
  for (int i = 0; i < G_k_w.size(); ++i) {
    EXPECT_NEAR(G_k_w_check(i).real(), G_k_w(i).real(), tol);
    EXPECT_NEAR(G_k_w_check(i).imag(), G_k_w(i).imag(), tol);
  }
}
