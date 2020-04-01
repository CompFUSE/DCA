// Copyright (C) 2018 UT-Battelle, LLC
// Copyright (C) 2018 ETH Zurich
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

template <typename Real>
using CtintWalkerTest =
    typename dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_INT>;

using namespace dca::phys::solver;
template <typename Real>
using Matrix = dca::linalg::Matrix<Real, dca::linalg::CPU>;
template <typename Real>
using MatrixPair = std::array<Matrix<Real>, 2>;

using dca::linalg::matrixop::determinantIP;
template <typename Real>
double computeDetRatio(MatrixPair<Real> a, MatrixPair<Real> b) {
  double res = 1;
  res *= determinantIP(a[0]) / determinantIP(b[0]);
  res *= determinantIP(a[1]) / determinantIP(b[1]);
  return res;
}

template <typename Real>
double determinant(MatrixPair<Real> a) {
  return determinantIP(a[0]) * determinantIP(a[1]);
}

using FloatingPointTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(CtintWalkerTest, FloatingPointTypes);

TYPED_TEST(CtintWalkerTest, InsertAndRemoveVertex) {
  using Real = TypeParam;

  // Setup
  std::vector<double> rng_values(1000);
  for (auto& x : rng_values)
    x = Real(std::rand()) / RAND_MAX;
  typename CtintWalkerTest<Real>::RngType rng(rng_values);

  auto& data = *CtintWalkerTest<Real>::data_;
  auto& parameters = CtintWalkerTest<Real>::parameters_;

  ctint::G0Interpolation<dca::linalg::CPU, Real> g0(
      dca::phys::solver::ctint::details::shrinkG0(data.G0_r_t));
  typename CtintWalkerTest<Real>::LabelDomain label_dmn;

  using Parameters = typename CtintWalkerTest<Real>::Parameters;
  using Walker = testing::phys::solver::ctint::WalkerWrapper<Parameters, Real>;
  Walker::setDMatrixBuilder(g0);
  Walker::setDMatrixAlpha(parameters.getAlphas(), 0);
  Walker::setInteractionVertices(data, parameters);

  Walker walker(parameters, rng);

  // *******************************
  // Test vertex removal ***********
  // *******************************
  // Set rng value to select: last vertex, unused, unused, accept
  rng.setNewValues(std::vector<double>{0.95, -1, -1, 0.01});
  MatrixPair<Real> old_M(walker.getM());
  bool result = walker.tryVertexRemoval();
  MatrixPair<Real> new_M(walker.getM());
  ASSERT_EQ(true, result);
  //  ASSERT_EQ(old_M.nrCols(), new_M.nrCols() + 2);
  // Compute directly the new M.
  walker.setMFromConfig();
  MatrixPair<Real> direct_M(walker.getM());
  for (int s = 0; s < 2; ++s)
    for (int j = 0; j < new_M[s].nrCols(); j++)
      for (int i = 0; i < new_M[s].nrRows(); i++)
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 10 * std::numeric_limits<Real>::epsilon());
  // Compute directly the determinant ratio. Note: M = D^-1.
  Real det_ratio = computeDetRatio(old_M, new_M);

  EXPECT_NEAR(det_ratio, walker.getRatio(), 100 * std::numeric_limits<Real>::epsilon());

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
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 10 * std::numeric_limits<Real>::epsilon());
  det_ratio = computeDetRatio(old_M, new_M);
  EXPECT_NEAR(det_ratio, walker.getRatio(), 100 * std::numeric_limits<Real>::epsilon());

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
        EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 10 * std::numeric_limits<Real>::epsilon());
}

TYPED_TEST(CtintWalkerTest, Inverse) {
  using Real = TypeParam;
  {
    Matrix<Real> M(std::make_pair(2, 2), std::make_pair(2, 2));
    M(0, 0) = 1, M(0, 1) = 2;
    M(1, 0) = 3, M(1, 1) = -1;

    Matrix<Real> inv(2, 3);
    const Real det = dca::phys::solver::ctint::details::smallDeterminant(M);
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<Real, dca::linalg::CPU> work;
    dca::phys::solver::ctint::details::smallInverse(M, inv, det, ipiv, work);
    constexpr Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
    EXPECT_NEAR(1. / 7., inv(0, 0), tolerance);
    EXPECT_NEAR(2. / 7., inv(0, 1), tolerance);
    EXPECT_NEAR(3. / 7., inv(1, 0), tolerance);
    EXPECT_NEAR(-1. / 7., inv(1, 1), tolerance);
  }
  // 4x4 case
  {
    Matrix<Real> M4(4);
    for (int j = 0; j < 4; ++j)
      for (int i = 0; i < 4; ++i)
        M4(i, j) = i + j * j;

    const Real det = dca::phys::solver::ctint::details::smallDeterminant(M4);
    EXPECT_NEAR(dca::linalg::matrixop::determinantIP(M4), det,
                10 * std::numeric_limits<Real>::epsilon());
  }
}

TEST(CtintWalkerTest, RemoveIndices) {
  std::vector<int> mock_configuration{0, 1, 2, 3, 4, 5, 6, 7};
  const int n = mock_configuration.size();

  Matrix<double> Q(std::make_pair(n, 2));
  Matrix<double> R(std::make_pair(2, n));

  // Fill the two rectangular matrices wth an arbitrary function of a 1D configuration.
  auto set_matrices_from_config = [](Matrix<double>& Q, Matrix<double>& R,
                                     const std::vector<int>& config) {
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < Q.nrRows(); ++i) {
        Q(i, j) = R(j, i) = config[i] + 0.5 * j;
      }
  };

  set_matrices_from_config(Q, R, mock_configuration);

  std::vector<unsigned short> remove{2, 6, 1, 1, 3};
  for (unsigned short index : remove) {
    std::swap(mock_configuration[index], mock_configuration.back());
    mock_configuration.pop_back();
  }

  Matrix<double> Q_expected(std::make_pair(mock_configuration.size(), 2));
  Matrix<double> R_expected(std::make_pair(2, mock_configuration.size()));
  set_matrices_from_config(Q_expected, R_expected, mock_configuration);

  dca::phys::solver::ctint::details::removeIndices(Q, R, remove);

  EXPECT_EQ(Q_expected, Q);
  EXPECT_EQ(R_expected, R);
}
