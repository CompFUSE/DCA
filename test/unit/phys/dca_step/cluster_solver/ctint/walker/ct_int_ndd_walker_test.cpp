// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using G0Setup = typename dca::testing::G0Setup<dca::testing::LatticeHund, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;
using Walker = testing::phys::solver::ctint::WalkerWrapper<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

double computeDetRatio(MatrixPair a, MatrixPair b);

TEST_F(G0Setup, RemoveAndInstertNddVertex) {
  // Setup
  parameters.setDoubleUpdateProb(1);
  std::vector<double> rng_values(1000);
  for (double& x : rng_values)
    x = double(std::rand()) / RAND_MAX;
  G0Setup::RngType rng(rng_values);

  ctint::G0Interpolation<dca::linalg::CPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  ctint::DMatrixBuilder<dca::linalg::CPU> builder(g0, RDmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters.getAlphas());
  Walker walker(parameters, rng, G0Setup::interaction_vertices, builder);

  // *******************************
  // Test vertex insertion *********
  // *******************************
  // Set rng values to select: last interaction vertex, tau, aux_spin,
  // double insertion, tau, aux_spin,  accept.
  rng.setNewValues(std::vector<double>{0., 0.4, 0, 0, 0.67, 0.51, -1});
  auto old_M = walker.getM();
  bool result = walker.tryVertexInsert();
  auto new_M = walker.getM();
  ASSERT_EQ(true, result);
  ASSERT_EQ(old_M[0].nrCols(), new_M[0].nrCols() - 2);
  // Compute directly the new M.
  walker.setMFromConfig();
  auto direct_M = walker.getM();

  for (int s = 0; s < 2; ++s)
    EXPECT_TRUE(dca::linalg::matrixop::areNear(new_M[s], direct_M[s], 1e-7));

  double det_ratio = computeDetRatio(old_M, new_M);
  EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);

  // *******************************
  // Test vertex removal ***********
  // *******************************
  // Set rng value to select: last vertex, double removal, first partner, accept
  rng.setNewValues(std::vector<double>{0.99, 0, 0, 0.01});
  walker.setMFromConfig();
  old_M = walker.getM();
  result = walker.tryVertexRemoval();
  new_M = walker.getM();
  ASSERT_EQ(true, result);
  ASSERT_EQ(old_M[0].nrCols(), new_M[0].nrCols() + 2);
  // Compute directly the new M.
  walker.setMFromConfig();
  direct_M = walker.getM();
  for (int s = 0; s < 2; ++s)
    EXPECT_TRUE(dca::linalg::matrixop::areNear(new_M[s], direct_M[s], 1e-7));

  // Compute directly the determinant ratio. Note: M = D^-1.
  det_ratio = computeDetRatio(old_M, new_M);

  EXPECT_LE(std::abs((walker.getRatio() - det_ratio) / det_ratio), 5e-7);
}

double computeDetRatio(MatrixPair a, MatrixPair b) {
  double res = 1;
  using dca::linalg::matrixop::determinantIP;
  res *= determinantIP(a[0]) / determinantIP(b[0]);
  res *= determinantIP(a[1]) / determinantIP(b[1]);
  return res;
}
