// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the exact diagonalization (ED) solver on a single-site Hubbard model and compares
// the free Green's functions, G_0(i\omega) and G_0(\tau), and the interacting Green's function,
// G(i\omega), against their analytic result.

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/ed_cluster_solver.hpp"

#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "ed_cluster_solver_test_helper.hpp"

TEST(EDClusterSolverSingleSiteTest, ComputeG0AndG) {
  dca::parallel::NoConcurrency concurrency(0, nullptr);

  dca::testing::Parameters params("", concurrency);
  params.read_input_and_broadcast<dca::io::JSONReader>(
      DCA_SOURCE_DIR "/test/integration/exact_diagonalization_advanced/single_site_input.json");
  params.update_model();
  params.update_domains();

  dca::testing::Data data(params);
  dca::testing::DataRealFreq data_real(params);
  data.initialize_H_0_and_H_i();

  dca::phys::solver::EDClusterSolver<dca::linalg::CPU, dca::testing::Parameters, dca::testing::Data>
      ed_solver(params, data, data_real);
  ed_solver.initialize(0);  // 0 = index of iteration.
  ed_solver.execute();

  // For the given system, the dispersion relation \epsilon_\vec{k} is given by the diagonal
  // elements of the non-interacting Hamiltonian.
  dca::func::function<double, dca::func::dmn_variadic<dca::testing::OrbitalSpinDmn, dca::testing::KDmn>> eps_k;
  for (int k = 0; k < dca::testing::KDmn::dmn_size(); ++k)
    for (int m = 0; m < dca::testing::OrbitalSpinDmn::dmn_size(); ++m)
      eps_k(m, k) = data.H_DCA(m, m, k).real();  // The elements are real.

  dca::testing::F_k_w G0_k_w_analytic;
  dca::testing::computeAnalyticG0_k_w(eps_k, params.get_chemical_potential(), G0_k_w_analytic);
  dca::testing::F_k_t G0_k_t_analytic;
  dca::testing::computeAnalyticG0_k_t(eps_k, params.get_chemical_potential(), params.get_beta(),
                                      G0_k_t_analytic);
  dca::testing::F_k_w G_w_analytic;
  dca::testing::computeAnalyticG_w(params.get_U(), G_w_analytic);

  const double tol = 1.e-14;
  for (int i = 0; i < G0_k_w_analytic.size(); ++i) {
    EXPECT_NEAR(G0_k_w_analytic(i).real(), data.G0_k_w(i).real(), tol);
    EXPECT_NEAR(G0_k_w_analytic(i).imag(), data.G0_k_w(i).imag(), tol);

    EXPECT_NEAR(G_w_analytic(i).real(), data.G_k_w(i).real(), tol);
    EXPECT_NEAR(G_w_analytic(i).imag(), data.G_k_w(i).imag(), tol);
  }
  for (int i = 0; i < G0_k_t_analytic.size(); ++i) {
    EXPECT_NEAR(G0_k_t_analytic(i), data.G0_k_t(i).real(), tol);
    EXPECT_NEAR(0, data.G0_k_t(i).imag(), tol);
  }
}
