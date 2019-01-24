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

using G0Setup =
    typename dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;
using Walker = testing::phys::solver::ctint::WalkerWrapper<G0Setup::Parameters>;
using Matrix = Walker::Matrix;
using MatrixPair = std::array<Matrix, 2>;

double computeDetRatio(MatrixPair a, MatrixPair b);
double determinant(MatrixPair a);

TEST_F(G0Setup, RemoveAndInstertVertex) {
  // Setup
  //  for(int i = 0; i < 900; ++i)
  //    std::rand();
  std::vector<double> rng_values(1000);
  for (double& x : rng_values)
    x = double(std::rand()) / RAND_MAX;
  G0Setup::RngType rng(rng_values);

  ctint::G0Interpolation<dca::linalg::CPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t));
  G0Setup::LabelDomain label_dmn;
  ctint::DMatrixBuilder<dca::linalg::CPU> builder(g0, RDmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters_.getAlphas());
  Walker walker(parameters_, rng, G0Setup::interaction_vertices_, builder);

  // *******************************
  // Test vertex removal ***********
  // *******************************
  // Set rng value to select: last vertex,  accept
  rng.setNewValues(std::vector<double>{0.95, 0.01});
  MatrixPair old_M(walker.getM());
  bool result = walker.tryVertexRemoval();
  MatrixPair new_M(walker.getM());
  ASSERT_EQ(true, result);
  //  ASSERT_EQ(old_M.nrCols(), new_M.nrCols() + 2);
  // Compute directly the new M.
  walker.setMFromConfig();
  MatrixPair direct_M(walker.getM());
  for (int s = 0; s < 2; ++s)
    for (int j = 0; j < new_M[s].nrCols(); j++)
      for (int i = 0; i < new_M[s].nrRows(); i++)
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 1e-7);
  // Compute directly the determinant ratio. Note: M = D^-1.
  double det_ratio = computeDetRatio(old_M, new_M);

  EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);

  // *******************************
  // Test vertex insertion *********
  // *******************************
  // Set rng values to select: first interaction vertex, tau, aux_spin,  accept
  rng.setNewValues(std::vector<double>{0, 0.4, 0.51, 1e-6});
  old_M = walker.getM();
  result = walker.tryVertexInsert();
  new_M = walker.getM();
  ASSERT_EQ(true, result);
  //  ASSERT_EQ(old_M.nrCols(), new_M.nrCols() - 2);
  // Compute directly the new M.
  walker.setMFromConfig();
  direct_M = walker.getM();
  for (int s = 0; s < 2; ++s)
    for (int j = 0; j < new_M[s].nrCols(); j++)
      for (int i = 0; i < new_M[s].nrRows(); i++)
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 1e-7);
  det_ratio = computeDetRatio(old_M, new_M);
  EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);

  // ****************************************
  // Test last vertex removal and insertion *
  // ****************************************
  rng.setNewValues(std::vector<double>(100, 0));
  const int n_vertices = walker.order();
  for (int i = 0; i < n_vertices - 1; i++)
    walker.tryVertexRemoval();
  ASSERT_EQ(1, walker.order());
  old_M = walker.getM();
  result = walker.tryVertexRemoval();
  // walker.getM() is now empty
  det_ratio = determinant(old_M) / 1.;
  EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);
  // Test insertion.
  rng.setNewValues(std::vector<double>{0, 0.5, 0, 1e-6});
  result = walker.tryVertexInsert();
  new_M = walker.getM();
  det_ratio = 1. / determinant(new_M);
  EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);
  walker.setMFromConfig();
  direct_M = walker.getM();
  for (int s = 0; s < 2; ++s)
    for (int j = 0; j < new_M[s].nrCols(); j++)
      for (int i = 0; i < new_M[s].nrRows(); i++)
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 1e-7);
}

TEST(WalkerMethodTest, Inverse) {
  {
    Matrix M(std::make_pair(2, 2), std::make_pair(2, 2));
    M(0, 0) = 1, M(0, 1) = 2;
    M(1, 0) = 3, M(1, 1) = -1;

    Matrix inv(2, 3);
    const double det = dca::phys::solver::ctint::details::smallDeterminant(M);
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<double, dca::linalg::CPU> work;
    dca::phys::solver::ctint::details::smallInverse(M, inv, det, ipiv, work);
    EXPECT_NEAR(1. / 7., inv(0, 0), 1e-12);
    EXPECT_NEAR(2. / 7., inv(0, 1), 1e-12);
    EXPECT_NEAR(3. / 7., inv(1, 0), 1e-12);
    EXPECT_NEAR(-1. / 7., inv(1, 1), 1e-12);
  }
  // 4x4 case
  {
    Matrix M4(4);
    for (int j = 0; j < 4; ++j)
      for (int i = 0; i < 4; ++i)
        M4(i, j) = i + j * j;

    const double det = dca::phys::solver::ctint::details::smallDeterminant(M4);
    EXPECT_DOUBLE_EQ(dca::linalg::matrixop::determinantIP(M4), det);
  }
}

TEST(WalkerMethodTest, RemoveIndices) {
  std::vector<int> mock_configuration{0, 1, 2, 3, 4, 5, 6, 7};
  const int n = mock_configuration.size();

  dca::linalg::Matrix<double, dca::linalg::CPU> Q(std::make_pair(n, 2));
  dca::linalg::Matrix<double, dca::linalg::CPU> R(std::make_pair(2, n));

  auto set_matrices_from_config = [](dca::linalg::Matrix<double, dca::linalg::CPU>& Q,
                                     dca::linalg::Matrix<double, dca::linalg::CPU>& R,
                                     const std::vector<int>& config) {
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < Q.nrRows(); ++i) {
        Q(i, j) = R(j, i) = config[i] + 0.5 * j;
      }
  };

  set_matrices_from_config(Q, R, mock_configuration);

  std::vector<ushort> remove{2, 6, 1, 1, 3};
  for (ushort index : remove) {
    std::swap(mock_configuration[index], mock_configuration.back());
    mock_configuration.pop_back();
  }

  dca::linalg::Matrix<double, dca::linalg::CPU> Q_expected(
      std::make_pair(mock_configuration.size(), 2));
  dca::linalg::Matrix<double, dca::linalg::CPU> R_expected(
      std::make_pair(2, mock_configuration.size()));
  set_matrices_from_config(Q_expected, R_expected, mock_configuration);

  dca::phys::solver::ctint::details::removeIndices(Q, R, remove);

  EXPECT_EQ(Q_expected, Q);
  EXPECT_EQ(R_expected, R);
}

using dca::linalg::matrixop::determinantIP;
double computeDetRatio(MatrixPair a, MatrixPair b) {
  double res = 1;
  res *= determinantIP(a[0]) / determinantIP(b[0]);
  res *= determinantIP(a[1]) / determinantIP(b[1]);
  return res;
}

double determinant(MatrixPair a) {
  return determinantIP(a[0]) * determinantIP(a[1]);
}
