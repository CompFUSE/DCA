// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper.hpp"
#include "walker_wrapper_submatrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/multinsert_input.json";

using G0Setup =
    typename dca::testing::G0Setup<dca::testing::LatticeHund, dca::phys::solver::CT_INT, input_name>;
using namespace dca::phys::solver;
using Walker = testing::phys::solver::ctint::WalkerWrapper<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

auto getVertexRng = [](int id) {
  const double n_vertices = 12.;
  return (n_vertices - id - 0.5) / n_vertices;
};

// Compare single versus double update without a submatrix update.
TEST_F(G0Setup, NoSubmatrix) {
  G0Setup::RngType rng(std::vector<double>{});

  ctint::G0Interpolation<dca::linalg::CPU, double> g0(
      dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  Walker::setDMatrixBuilder(g0);
  Walker::setDMatrixAlpha(parameters_.getAlphas(), false);
  Walker::setInteractionVertices(*data_, parameters_);

  parameters_.setDoubleUpdateProbability(0);
  Walker walker_single(parameters_, rng);

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0.41, 0, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(8), 0.53, 1, 0});
  walker_single.tryVertexInsert();

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(7), 0.34, 0, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(9), 0.36, 0, 0});
  walker_single.tryVertexInsert();

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0.24, 1, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(8), 0.19, 1, 0});
  walker_single.tryVertexInsert();

  // interaction, accept
  rng.setNewValues(std::vector<double>{0, 0});
  walker_single.tryVertexRemoval();
  rng.setNewValues(std::vector<double>{0, 0});
  walker_single.tryVertexRemoval();

  auto M1 = walker_single.getM();
  walker_single.setMFromConfig();
  auto M1_dir = walker_single.getM();
  constexpr int final_size = 4;
  ASSERT_EQ(M1[0].nrCols(), final_size);
  ASSERT_EQ(M1[1].nrCols(), final_size);

  parameters_.setDoubleUpdateProbability(1);
  Walker walker_double(parameters_, rng);

  ////////////////////////////////

  // interaction, partner, tau, aux, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0});
  walker_double.tryVertexInsert();
  // interaction, partner tau, aux, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(7), 0.9, 0.34, 0, 0.36, 0, 0});
  walker_double.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.24, 1, 0.19, 1, 0});
  walker_double.tryVertexInsert();
  //
  // first_id, double removal, partner_id, accept
  rng.setNewValues(std::vector<double>{0, 0., 0, 0, 0});
  walker_double.tryVertexRemoval();

  auto M2 = walker_double.getM();
  walker_double.setMFromConfig();
  auto M2_dir = walker_double.getM();
  ASSERT_EQ(M2[0].nrCols(), final_size);
  ASSERT_EQ(M2[1].nrCols(), final_size);

  for (int s = 0; s < 2; ++s) {
    EXPECT_TRUE(dca::linalg::matrixop::areNear(M1[s], M1_dir[s], 1e-7));
    EXPECT_TRUE(dca::linalg::matrixop::areNear(M2[s], M2_dir[s], 1e-7));

    ASSERT_EQ(M1[s].size(), M2[s].size());
    std::vector<double> el1, el2;
    for (int j = 0; j < M1[s].nrRows(); ++j)
      for (int i = 0; i < M1[s].nrCols(); ++i) {
        el1.push_back(M1[s](i, j));
        el2.push_back(M2[s](i, j));
      }
    std::sort(el1.begin(), el1.end());
    std::sort(el2.begin(), el2.end());

    for (int i = 0; i < el1.size(); ++i)
      EXPECT_NEAR(el1[i], el2[i], 1e-5);
  }
}

// Compare the double insertion probability of a submatrix walker with a non-submatrix one.
using WalkerSubmatrix = testing::phys::solver::ctint::WalkerWrapperSubmatrix<G0Setup::Parameters>;
TEST_F(G0Setup, Submatrix) {
  G0Setup::RngType rng(std::vector<double>{});
  ctint::G0Interpolation<dca::linalg::CPU, double> g0(
      dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  Walker::setDMatrixBuilder(g0);
  Walker::setDMatrixAlpha(parameters_.getAlphas(), false);
  Walker::setInteractionVertices(*data_, parameters_);
  parameters_.setDoubleUpdateProbability(1);

  Walker walker(parameters_, rng);

  // interaction, partner, tau, aux, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0});
  walker.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(7), 0.9, 0.34, 0, 0.36, 0, 0});
  walker.tryVertexInsert();
  // first idx, double rem, second partner, accept
  rng.setNewValues(std::vector<double>{0, 0, 0, 0});
  walker.tryVertexRemoval();

  auto M1 = walker.getM();
  const double prob1 = walker.getAcceptanceProbability();

  /////////////////////////////

  WalkerSubmatrix::setDMatrixBuilder(g0);
  WalkerSubmatrix::setDMatrixAlpha(parameters_.getAlphas(), false);
  WalkerSubmatrix::setInteractionVertices(*data_, parameters_);

  WalkerSubmatrix walker_subm(parameters_, rng);

  std::vector<double> random_vals{// (insert, first_id, partner_id, tau, aux, tau, aux, accept) x 2
                                  0, getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0, 0, getVertexRng(7),
                                  0.9, 0.34, 0, 0.36, 0, 0,
                                  // remove, first_id, double_removal, second_id, accept
                                  1, 0, 0, 0, 0};

  rng.setNewValues(random_vals);
  walker_subm.doStep(3);

  auto M2 = walker_subm.getM();
  const double prob2 = walker_subm.getAcceptanceProbability();

  EXPECT_NEAR(prob1, prob2, 1e-7);
  for (int s = 0; s < 2; ++s)
    EXPECT_TRUE(dca::linalg::matrixop::areNear(M1[s], M2[s], 1e-7));
}
