// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_free_greens_function.hpp.

#include "dca/phys/dca_algorithms/compute_free_greens_function.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

using namespace dca;

class ComputeFreeGreensFunctionTest : public ::testing::Test {
protected:
  using PointGroup = phys::domains::D4;
  using SpinDmn = func::dmn<2, int>;
  using MatsubaraFreqDmn = func::dmn<4, double>;  // Matsubara frequency domain with 4 elements.
  using ImagTimeDmn = func::dmn<5, double>;       // Imaginary time domain with 5 elements.

  ComputeFreeGreensFunctionTest() : concurrency_(0, nullptr) {}

  static void SetUpTestCase() {
    std::vector<double> freqs(4);
    freqs[0] = -3. * M_PI / beta;
    freqs[1] = -M_PI / beta;
    freqs[2] = M_PI / beta;
    freqs[3] = 3. * M_PI / beta;
    MatsubaraFreqDmn::set_elements(freqs);

    std::vector<double> t_elements(5);
    t_elements[0] = -beta;
    t_elements[1] = -beta / 2.;
    t_elements[2] = 0;
    t_elements[3] = beta / 2.;
    t_elements[4] = beta;
    ImagTimeDmn::set_elements(t_elements);
  }

  static constexpr double beta = 1.;  // inverse temperature

  const parallel::NoConcurrency concurrency_;
};

// Test for a 2x2 square lattice with diagonal (in orbital-spin space) non-interacting Hamiltonian
// H_0.
TEST_F(ComputeFreeGreensFunctionTest, SquareLattice) {
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
      H_0;
  Lattice::initialize_H_0(params, H_0);

  const double mu = 0.;  // chemical potential

  //
  // G_0(\vec{k}, i\omega_n)
  //
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G0_k_w;

  phys::compute_G0_k_w(H_0, mu, 3, G0_k_w);

  // Check spin symmetry and that off-diagonal elements vanish.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, k, wn).real(), G0_k_w(1, 1, k, wn).real());
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, k, wn).imag(), G0_k_w(1, 1, k, wn).imag());

      EXPECT_DOUBLE_EQ(0., G0_k_w(0, 1, k, wn).real());
      EXPECT_DOUBLE_EQ(0., G0_k_w(0, 1, k, wn).imag());
      EXPECT_DOUBLE_EQ(0., G0_k_w(1, 0, k, wn).real());
      EXPECT_DOUBLE_EQ(0., G0_k_w(1, 0, k, wn).imag());
    }
  }

  // Check symmetry in Matsubara frequency.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size() / 2; ++wn) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, k, wn).real(),
                       G0_k_w(0, 0, k, MatsubaraFreqDmn::get_size() - 1 - wn).real());
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, k, wn).imag(),
                       -G0_k_w(0, 0, k, MatsubaraFreqDmn::get_size() - 1 - wn).imag());
    }
  }

  // Check momentum space symmetry (0, \pi) <--> (\pi, 0).
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 1, wn).real(), G0_k_w(0, 0, 2, wn).real());
    EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 1, wn).imag(), G0_k_w(0, 0, 2, wn).imag());
  }

  // Check momentum space symmetry (0, 0) <--> (\pi, \pi).
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 0, wn).real(), -G0_k_w(0, 0, 3, wn).real());
    EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 0, wn).imag(), G0_k_w(0, 0, 3, wn).imag());
  }

  // Check that Re[G_0(\vec{k}={0, \pi}] = 0.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    EXPECT_DOUBLE_EQ(0., G0_k_w(0, 0, 1, wn).real());
  }

  // Check some system specific values of G0.
  EXPECT_DOUBLE_EQ(0.15462161453970907, G0_k_w(0, 0, 0, 1).real());
  EXPECT_DOUBLE_EQ(0.12143953208103569, G0_k_w(0, 0, 0, 1).imag());
  EXPECT_DOUBLE_EQ(0., G0_k_w(0, 0, 1, 0).real());
  EXPECT_DOUBLE_EQ(0.1061032953945969, G0_k_w(0, 0, 1, 0).imag());
  EXPECT_DOUBLE_EQ(0., G0_k_w(0, 0, 2, 2).real());
  EXPECT_DOUBLE_EQ(-0.31830988618379069, G0_k_w(0, 0, 2, 2).imag());
  EXPECT_DOUBLE_EQ(-0.038158312109895294, G0_k_w(0, 0, 3, 3).real());
  EXPECT_DOUBLE_EQ(-0.089908404748375109, G0_k_w(0, 0, 3, 3).imag());

  //
  // G_0(\vec{k}, \tau)
  //
  func::function<double, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                            func::dmn_0<ImagTimeDmn>>>
      G0_k_t;

  phys::compute_G0_k_t(H_0, mu, beta, G0_k_t);

  // Check spin symmetry and that off-diagonal elements vanish.
  for (int t = 0; t < ImagTimeDmn::get_size(); ++t) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_k_t(0, 0, k, t), G0_k_t(1, 1, k, t));

      EXPECT_DOUBLE_EQ(0., G0_k_t(0, 1, k, t));
      EXPECT_DOUBLE_EQ(0., G0_k_t(1, 0, k, t));
    }
  }

  // Check symmetry in imaginary time: G_0(\tau+\beta) = -G_0(\tau).
  for (int k = 0; k < KDmn::get_size(); ++k) {
    // -beta --> 0
    EXPECT_DOUBLE_EQ(G0_k_t(0, 0, k, 0), -G0_k_t(0, 0, k, 2));

    // -beta/2 --> beta/2
    EXPECT_DOUBLE_EQ(G0_k_t(0, 0, k, 1), -G0_k_t(0, 0, k, 3));
  }

  // Check momentum space symmetry (0, \pi) <--> (\pi, 0).
  for (int t = 0; t < ImagTimeDmn::get_size(); ++t) {
    EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 1, t), G0_k_t(0, 0, 2, t));
  }

  // Check that G_0(\vec{0, \pi}, \tau>=0) = const.
  EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 1, 2), G0_k_t(0, 0, 1, 3));
  EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 1, 2), G0_k_t(0, 0, 1, 4));

  // Check that G_0(\vec{0, \pi}, \tau<0) = const.
  EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 1, 0), G0_k_t(0, 0, 1, 1));

  // Check some system specific values of G0.
  EXPECT_DOUBLE_EQ(0.13290111441703986, G0_k_t(0, 0, 0, 1));
  EXPECT_DOUBLE_EQ(0.5, G0_k_t(0, 0, 1, 0));
  EXPECT_DOUBLE_EQ(-0.5, G0_k_t(0, 0, 2, 2));
  EXPECT_DOUBLE_EQ(-0.98201379003790845, G0_k_t(0, 0, 3, 2));
}

// Test for a single-site bilayer lattice with off-diagonal elements (in orbital-spin space) in the
// non-interacting Hamiltonian H_0.
TEST_F(ComputeFreeGreensFunctionTest, BilayerLattice) {
  using Lattice = phys::models::bilayer_lattice<PointGroup>;
  using OrbitalDmn = func::dmn<2, int>;  // 2 orbitals
  using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<OrbitalDmn>, func::dmn_0<SpinDmn>>;

  // Momentum space domain of the single-site bilayer lattice
  using KDmn = func::dmn<1, std::vector<double>>;
  const std::vector<std::vector<double>> k_vecs{{0., 0.}};
  KDmn::set_elements(k_vecs);

  phys::params::ModelParameters<phys::models::TightBindingModel<Lattice>> params;
  params.set_t_perp(1.);

  func::function<std::complex<double>,
                 func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>>>
      H_0;
  Lattice::initialize_H_0(params, H_0);

  const double mu = 0.9;  // chemical potential

  //
  // G_0(\vec{k}, i\omega_n)
  //
  func::function<std::complex<double>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                                          func::dmn_0<MatsubaraFreqDmn>>>
      G0_k_w;

  phys::compute_G0_k_w(H_0, mu, 3, G0_k_w);

  // Check spin symmetry and that off-diagonal (in spin) elements vanish.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      for (int o2 = 0; o2 < OrbitalDmn::get_size(); ++o2) {
        for (int o1 = 0; o1 < OrbitalDmn::get_size(); ++o1) {
          EXPECT_DOUBLE_EQ(G0_k_w(o1, 0, o2, 0, k, wn).real(), G0_k_w(o1, 1, o2, 1, k, wn).real());
          EXPECT_DOUBLE_EQ(G0_k_w(o1, 0, o2, 0, k, wn).imag(), G0_k_w(o1, 1, o2, 1, k, wn).imag());

          EXPECT_DOUBLE_EQ(0., G0_k_w(o1, 0, o2, 1, k, wn).real());
          EXPECT_DOUBLE_EQ(0., G0_k_w(o1, 0, o2, 1, k, wn).imag());
          EXPECT_DOUBLE_EQ(0., G0_k_w(o1, 1, o2, 0, k, wn).real());
          EXPECT_DOUBLE_EQ(0., G0_k_w(o1, 1, o2, 0, k, wn).imag());
        }
      }
    }
  }

  // Check orbital symmetry.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size(); ++wn) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 0, 0, k, wn).real(), G0_k_w(1, 0, 1, 0, k, wn).real());
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 0, 0, k, wn).imag(), G0_k_w(1, 0, 1, 0, k, wn).imag());

      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 1, 0, k, wn).real(), G0_k_w(1, 0, 0, 0, k, wn).real());
      EXPECT_DOUBLE_EQ(G0_k_w(0, 0, 1, 0, k, wn).imag(), G0_k_w(1, 0, 0, 0, k, wn).imag());
    }
  }

  // Check symmetry in Matsubara frequency.
  for (int wn = 0; wn < MatsubaraFreqDmn::get_size() / 2; ++wn) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      for (int o2 = 0; o2 < OrbitalDmn::get_size(); ++o2) {
        for (int o1 = 0; o1 < OrbitalDmn::get_size(); ++o1) {
          EXPECT_DOUBLE_EQ(G0_k_w(o1, 0, o2, 0, k, wn).real(),
                           G0_k_w(o1, 0, o2, 0, k, MatsubaraFreqDmn::get_size() - 1 - wn).real());
          EXPECT_DOUBLE_EQ(G0_k_w(o1, 0, o2, 0, k, wn).imag(),
                           -G0_k_w(o1, 0, o2, 0, k, MatsubaraFreqDmn::get_size() - 1 - wn).imag());
        }
      }
    }
  }

  // Check some system specific values of G0.
  EXPECT_DOUBLE_EQ(0.009714500128141767, G0_k_w(0, 0, 0, 0, 0, 0).real());
  EXPECT_DOUBLE_EQ(0.104025451807396321, G0_k_w(0, 0, 0, 0, 0, 0).imag());
  EXPECT_DOUBLE_EQ(0.065415914352591892, G0_k_w(0, 0, 0, 0, 0, 2).real());
  EXPECT_DOUBLE_EQ(-0.275525185917382232, G0_k_w(0, 0, 0, 0, 0, 2).imag());

  EXPECT_DOUBLE_EQ(0.075537777125657318, G0_k_w(0, 0, 1, 0, 0, 1).real());
  EXPECT_DOUBLE_EQ(-0.042462511367682185, G0_k_w(0, 0, 1, 0, 0, 1).imag());
  EXPECT_DOUBLE_EQ(0.010840164331246138, G0_k_w(0, 0, 1, 0, 0, 3).real());
  EXPECT_DOUBLE_EQ(0.0020658999190548404, G0_k_w(0, 0, 1, 0, 0, 3).imag());

  //
  // G_0(\vec{k}, \tau)
  //
  func::function<double, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, func::dmn_0<KDmn>,
                                            func::dmn_0<ImagTimeDmn>>>
      G0_k_t;

  phys::compute_G0_k_t(H_0, mu, beta, G0_k_t);

  // Check spin symmetry and that off-diagonal (in spin) elements vanish.
  for (int t = 0; t < ImagTimeDmn::get_size(); ++t) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      for (int o2 = 0; o2 < OrbitalDmn::get_size(); ++o2) {
        for (int o1 = 0; o1 < OrbitalDmn::get_size(); ++o1) {
          EXPECT_DOUBLE_EQ(G0_k_t(o1, 0, o2, 0, k, t), G0_k_t(o1, 1, o2, 1, k, t));

          EXPECT_DOUBLE_EQ(0., G0_k_t(o1, 0, o2, 1, k, t));
          EXPECT_DOUBLE_EQ(0., G0_k_t(o1, 1, o2, 0, k, t));
        }
      }
    }
  }

  // Check orbital symmetry.
  for (int t = 0; t < ImagTimeDmn::get_size(); ++t) {
    for (int k = 0; k < KDmn::get_size(); ++k) {
      EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 0, 0, k, t), G0_k_t(1, 0, 1, 0, k, t));
      EXPECT_DOUBLE_EQ(G0_k_t(0, 0, 1, 0, k, t), G0_k_t(1, 0, 0, 0, k, t));
    }
  }

  // Check symmetry in imaginary time: G_0(\tau+\beta) = -G_0(\tau).
  for (int k = 0; k < KDmn::get_size(); ++k) {
    for (int o2 = 0; o2 < OrbitalDmn::get_size(); ++o2) {
      for (int o1 = 0; o1 < OrbitalDmn::get_size(); ++o1) {
        // -beta --> 0
        EXPECT_DOUBLE_EQ(G0_k_t(o1, 0, o2, 0, k, 0), -G0_k_t(o1, 0, o2, 0, k, 2));

        // -beta/2 --> beta/2
        EXPECT_DOUBLE_EQ(G0_k_t(o1, 0, o2, 0, k, 1), -G0_k_t(o1, 0, o2, 0, k, 3));
      }
    }
  }

  // Check some system specific values of G0.
  EXPECT_DOUBLE_EQ(0.32754383092096884, G0_k_t(0, 0, 0, 0, 0, 0));
  EXPECT_DOUBLE_EQ(-0.417899194649848038, G0_k_t(0, 0, 0, 0, 0, 3));

  EXPECT_DOUBLE_EQ(-0.081476455730596378, G0_k_t(0, 0, 1, 0, 0, 1));
  EXPECT_DOUBLE_EQ(-0.19743535655797104, G0_k_t(0, 0, 1, 0, 0, 4));
}
