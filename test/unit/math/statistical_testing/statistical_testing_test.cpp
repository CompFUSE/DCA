// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Andrei Plamada (plamada@phys.ethz.ch)
//
// This file tests the statistical test framework provided by statistical_testing.hpp.

#include "dca/math/statistical_testing/statistical_testing.hpp"

#include "gtest/gtest.h"

#include "dca/io/csv/csv_reader.hpp"
#include "dca/function/domains.hpp"

using namespace dca::func;

template <class Domain>
void writeLT(function<double, dmn_variadic<Domain, Domain>>& f);
constexpr double tolerance = 1e-12;
const std::string test_directory = DCA_SOURCE_DIR "/test/unit/math/statistical_testing/";

using dca::math::StatisticalTesting;
using namespace dca;

TEST(ProbabilityDistributions, Chi2Cdf) {
  io::CSVReader reader;
  std::vector<std::vector<double>> data;
  // format of the reference data originated with python:
  // d.o.f.
  // x1 ... xn
  // y1 ... yn
  reader.execute(test_directory + "chi2_distribution_reference.csv", data);
  for (int i = 0; i < data.size(); i += 3) {
    const int dof = data[i][0];
    ASSERT_EQ(data[i + 1].size(), data[i + 2].size());
    for (int j = 0; j < data[i + 1].size(); ++j) {
      const double x = data[i + 1][j];
      const double y = data[i + 2][j];

      EXPECT_NEAR(y, dca::math::chi2Cdf(x, dof), 1e-8);
    }
  }
  // Test for invalid arguments.
  EXPECT_THROW(dca::math::chi2Cdf(-1, 1), std::logic_error);
  EXPECT_THROW(dca::math::chi2Cdf(1, 0), std::logic_error);
}

TEST(ProbabilityDistributions, FCdf) {
  io::CSVReader reader;
  std::vector<std::vector<double>> data;

  reader.execute(test_directory + "f_distribution_reference.csv", data);
  for (int i = 0; i < data.size(); i += 3) {
    const int d1 = data[i][0];
    const int d2 = data[i][1];
    ASSERT_EQ(data[i + 1].size(), data[i + 2].size());
    for (int j = 0; j < data[i + 1].size(); ++j) {
      const double x = data[i + 1][j];
      const double y = data[i + 2][j];

      EXPECT_NEAR(y, dca::math::fCdf(x, d1, d2), 1e-8);
    }
  }

  EXPECT_THROW(dca::math::fCdf(-1, 1, 1), std::logic_error);
  EXPECT_THROW(dca::math::fCdf(1, 1, 0), std::logic_error);
}

TEST(StatisticalTesting, ComputePValue1) {
  // Non singular covariance matrix.
  using Domain = dmn_0<dmn<3>>;
  function<double, Domain> f(""), f0("");
  f0 = 0;

  f(0) = 0.1;
  f(1) = 0.15;
  f(2) = 0.2;

  function<double, dmn_variadic<Domain, Domain>> cov("");
  cov(0, 0) = 0.1, cov(0, 1) = 0.1, cov(0, 2) = 0.1;
  cov(1, 1) = 0.2, cov(1, 2) = 0.05;
  cov(2, 2) = 0.3;
  writeLT(cov);

  StatisticalTesting test1(f, f0, cov);
  StatisticalTesting test2(test1);
  // Result obtained with a straightforward computation in Python.
  const double d2 = 0.21428571428571425;
  const int n = 5;
  const double expected_pval = 1. - dca::math::fCdf((n - 3.) / 3. * d2, 3, n - 3);
  EXPECT_NEAR(expected_pval, test1.computePValue(false, n, false), tolerance);
  EXPECT_EQ(3, test1.get_dof());
  // Performs the test with a faster implementation.
  EXPECT_NEAR(expected_pval, test2.computePValue(false, n, true), tolerance);
  EXPECT_EQ(3, test2.get_dof());
}

TEST(StatisticalTesting, ComputePvalue2) {
  // Singular covariance matrix.
  using Domain = dmn_0<dmn<4, int>>;
  function<double, Domain> f(""), f0("");
  f0 = 0;
  f(0) = 0.1;
  f(1) = 0.1;
  f(2) = 0.15;
  f(3) = 0.2;

  function<double, dmn_variadic<Domain, Domain>> cov("");
  cov(0, 0) = 0.1, cov(0, 1) = 0.1, cov(0, 2) = 0.1, cov(0, 3) = 0.1;
  cov(1, 1) = 0.1, cov(1, 2) = 0.1, cov(1, 3) = 0.1;
  cov(2, 2) = 0.2, cov(2, 3) = 0.05;
  cov(3, 3) = 0.3;
  writeLT(cov);

  // Once the first row and column is removed, the Mahalanobis distance is the same as the previous
  // test.
  const double d2 = 0.21428571428571425;
  const int n = 5;
  const double expected_pval = 1. - dca::math::chi2Cdf(n * d2, 3);
  StatisticalTesting test1(f, f0, cov);
  StatisticalTesting test2(test1);

  EXPECT_NEAR(expected_pval, test1.computePValue(true, n, false), tolerance);
  // As the covariance is singular setting the fast computation reverts to the default one.
  EXPECT_NEAR(expected_pval, test2.computePValue(true, n, true), tolerance);

  EXPECT_EQ(3, test1.get_dof());
  EXPECT_EQ(3, test2.get_dof());
}

TEST(StatisticalTesting, SelectIndices1) {
  using Domain = dmn_0<dmn<3>>;
  function<double, Domain> f(""), f0("");
  function<double, dmn_variadic<Domain, Domain>> cov("");

  f(0) = f(1) = 1;
  f(2) = 50.;
  cov(0, 0) = cov(1, 1) = cov(2, 2) = 1;
  cov(0, 1) = cov(1, 0) = 0.1;

  std::vector<int> indices{1, 0};
  StatisticalTesting test(f, f0, cov, 0);
  EXPECT_EQ(3, test.get_dof());
  // Test only the first two indices.
  test.selectIndices(indices);

  const double d2 = 2. * 0.9 / 0.99;
  const double expected = 1. - dca::math::chi2Cdf(d2, indices.size());
  EXPECT_NEAR(expected, test.computePValue(true, 1), 1e-10);
  const std::vector<int> new_indices{0, 1};
  EXPECT_EQ(new_indices, indices);
  EXPECT_EQ(2, test.get_dof());

  // Test if discardIndices obtains the same result.
  StatisticalTesting test2(f, f0, cov, 0);
  std::vector<int> discard{2};
  test2.discardIndices(discard);
  EXPECT_NEAR(expected, test2.computePValue(true, 1), 1e-10);
}

TEST(StatisticalTesting, SelectIndices2) {
  // Test exceptions.
  using Domain = dmn_0<dmn<3>>;
  function<double, Domain> f(""), f0("");
  function<double, dmn_variadic<Domain, Domain>> cov("");

  std::vector<int> illegal_indices{0, -1, 2};
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.selectIndices(illegal_indices), std::out_of_range);
  }
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.discardIndices(illegal_indices), std::out_of_range);
  }

  std::vector<int> illegal_indices2{5, 1, 2};
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.selectIndices(illegal_indices2), std::out_of_range);
  }
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.discardIndices(illegal_indices2), std::out_of_range);
  }

  {
    StatisticalTesting test(f, f0, cov);
    std::vector<int> redundant_indices{0, 1, 0, 2};
    EXPECT_NO_THROW(test.selectIndices(redundant_indices));
    const std::vector<int> expected{0, 1, 2};
    EXPECT_EQ(expected, redundant_indices);
    EXPECT_EQ(3, test.get_dof());
  }

  // empty tests are not allowed.
  std::vector<int> all_indices{0, 1, 2};
  std::vector<int> no_index{};
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.discardIndices(all_indices), std::logic_error);
  }
  {
    StatisticalTesting test(f, f0, cov);
    EXPECT_THROW(test.selectIndices(no_index), std::logic_error);
  }
}

TEST(StatisticalTesting, BigDifference) {
  // If the difference between f and f0 is large the pvalue must be close to zero.
  using Domain = dmn_0<dmn<2>>;
  function<double, Domain> f_wrong(""), f0("");
  f_wrong(0) = 0.5;
  f_wrong(1) = -0.1;
  f0(0) = 0.1;
  f0(1) = 0.1;

  function<double, dmn_variadic<Domain, Domain>> cov("");
  cov(0) = 0.2, cov(1) = 0.05;
  cov(2) = 0.05, cov(3) = 0.1;

  const int n_samples = 10;

  StatisticalTesting stat_test(f_wrong, f0, cov);
  StatisticalTesting stat_test_fast(f_wrong, f0, cov);
  EXPECT_GE(0.05, stat_test.computePValue(true, n_samples));

  EXPECT_NEAR(stat_test.computePValue(true, n_samples),
              stat_test_fast.computePValue(true, n_samples, true), 1e-10);
}

template <class Domain>
void writeLT(function<double, dmn_variadic<Domain, Domain>>& f) {
  for (int i = 0; i < Domain::dmn_size(); ++i)
    for (int j = i + 1; j < Domain::dmn_size(); ++j)
      f(j, i) = f(i, j);
}
