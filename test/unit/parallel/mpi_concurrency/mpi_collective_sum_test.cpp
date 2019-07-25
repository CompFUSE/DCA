// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Andrei Plamada (plamada@phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests mpi_collective_sum.hpp.
//
// This test only passes for 8 MPI processes.

#include "dca/parallel/mpi_concurrency/mpi_collective_sum.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/testing/minimalist_printer.hpp"

class MPICollectiveSumTest : public ::testing::Test {
protected:
  MPICollectiveSumTest() {
    size_ = sum_interface_.get_size();
    rank_ = sum_interface_.get_id();
  }

  int size_;
  int rank_;
  constexpr static double epsilon_ = 1e-14;

  dca::parallel::MPICollectiveSum sum_interface_;
};

TEST_F(MPICollectiveSumTest, SumScalar) {
  int scalar_test;
  int scalar_expected;

  scalar_test = rank_;
  sum_interface_.sum(scalar_test);

  scalar_expected = size_ * (size_ - 1) / 2;

  EXPECT_EQ(scalar_expected, scalar_test);
}

TEST_F(MPICollectiveSumTest, SumFunction) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;

  dca::func::function<double, TestDomain> function_test("test");
  dca::func::function<double, TestDomain> function_expected("expected");

  for (int i = 0; i < function_test.size(); i++)
    function_test(i) = i * rank_;

  sum_interface_.sum(function_test);

  for (int i = 0; i < function_test.size(); i++)
    function_expected(i) = i * size_ * (size_ - 1) / 2;

  for (int i = 0; i < function_test.size(); i++)
    EXPECT_EQ(function_expected(i), function_test(i));
}

TEST_F(MPICollectiveSumTest, LeaveOneOutAvgAndSum) {
  std::vector<double> values(size_);
  double sum = 0.;
  for (int i = 0; i < size_; ++i) {
    values[i] = 3.14 + i;
    sum += values[i];
  }

  // Expected result
  const double expected = (sum - values[rank_]) / double(size_ - 1);

  // Check scalar version.
  double scalar = values[rank_];
  double sum_one_out = scalar;
  sum_interface_.leaveOneOutAvg(scalar);
  sum_interface_.leaveOneOutSum(sum_one_out);
  EXPECT_DOUBLE_EQ(expected, scalar);
  EXPECT_DOUBLE_EQ(expected * (size_ - 1), sum_one_out);

  // Check dca::func::function version.
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  dca::func::function<double, TestDomain> f;
  f(0) = values[rank_];
  f(1) = 0.;

  dca::func::function<double, TestDomain> f_sum(f);
  sum_interface_.leaveOneOutAvg(f);
  sum_interface_.leaveOneOutSum(f_sum);

  EXPECT_DOUBLE_EQ(expected, f(0));
  EXPECT_DOUBLE_EQ(expected * (size_ - 1), f_sum(0));

  EXPECT_DOUBLE_EQ(0., f(1));
  EXPECT_DOUBLE_EQ(0., f_sum(1));
}

TEST_F(MPICollectiveSumTest, JackknifeErrorReal) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using FunctionType = dca::func::function<double, TestDomain>;

  FunctionType f;

  // Trivial case
  // All jackknife estimates are identical: jackknife error = 0
  f(0) = 3.14;
  f(1) = 2.72;

  auto err_trivial = sum_interface_.jackknifeError(f);

  EXPECT_NEAR(0., err_trivial(0), epsilon_);
  EXPECT_NEAR(0., err_trivial(1), epsilon_);

  // Non-trivial case
  const double d = 42.1;
  f(0) = rank_;
  f(1) = d;

  const FunctionType f_copy(f);

  FunctionType err_expected;
  err_expected(0) = std::sqrt(double(size_ - 1) / double(size_) * 2 * 21);
  err_expected(1) = 0.;

  // Do not overwrite the jackknife estimates with their average.
  auto err_no_overwriting = sum_interface_.jackknifeError(f, false);

  EXPECT_DOUBLE_EQ(f_copy(0), f(0));
  EXPECT_DOUBLE_EQ(f_copy(1), f(1));

  EXPECT_DOUBLE_EQ(err_expected(0), err_no_overwriting(0));
  EXPECT_DOUBLE_EQ(err_expected(1), err_no_overwriting(1));

  // Overwrite the jackknife estimates with their average.
  auto err_overwriting = sum_interface_.jackknifeError(f, true);

  const double rank_avg = double(size_ - 1) / 2.;

  EXPECT_DOUBLE_EQ(rank_avg, f(0));
  EXPECT_DOUBLE_EQ(d, f(1));

  EXPECT_DOUBLE_EQ(err_expected(0), err_overwriting(0));
  EXPECT_DOUBLE_EQ(err_expected(1), err_overwriting(1));
}

TEST_F(MPICollectiveSumTest, JackknifeErrorComplex) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using FunctionType = dca::func::function<std::complex<double>, TestDomain>;

  FunctionType f;

  // Trivial case
  // All jackknife estimates are identical: jackknife error = 0
  f(0) = std::complex<double>(3.14, 1.2);
  f(1) = std::complex<double>(2.72, 3.4);

  auto err_trivial = sum_interface_.jackknifeError(f);
  EXPECT_NEAR(0., err_trivial(0).real(), epsilon_);
  EXPECT_NEAR(0., err_trivial(1).real(), epsilon_);
  EXPECT_NEAR(0., err_trivial(0).imag(), epsilon_);
  EXPECT_NEAR(0., err_trivial(1).imag(), epsilon_);

  // Non-trivial case
  const double d = 42.1;
  const double r = 1.4;
  f(0) = std::complex<double>(rank_, rank_ + r);
  f(1) = std::complex<double>(rank_, d);

  const FunctionType f_copy(f);

  FunctionType err_expected;
  const double err_tmp = std::sqrt(double(size_ - 1) / double(size_) * 2 * 21);
  err_expected(0) = std::complex<double>(err_tmp, err_tmp);
  err_expected(1) = std::complex<double>(err_tmp, 0.);

  // Do not overwrite the jackknife estimates with their average.
  auto err_no_overwriting = sum_interface_.jackknifeError(f, false);

  EXPECT_DOUBLE_EQ(f_copy(0).real(), f(0).real());
  EXPECT_DOUBLE_EQ(f_copy(1).real(), f(1).real());

  EXPECT_DOUBLE_EQ(err_expected(0).real(), err_no_overwriting(0).real());
  EXPECT_DOUBLE_EQ(err_expected(0).imag(), err_no_overwriting(0).imag());
  EXPECT_DOUBLE_EQ(err_expected(1).real(), err_no_overwriting(1).real());
  EXPECT_DOUBLE_EQ(err_expected(1).imag(), err_no_overwriting(1).imag());

  // Overwrite the jackknife estimates with their average.
  auto err_overwriting = sum_interface_.jackknifeError(f, true);

  const double rank_avg = double(size_ - 1) / 2.;

  EXPECT_DOUBLE_EQ(std::complex<double>(rank_avg, rank_avg + r).real(), f(0).real());
  EXPECT_DOUBLE_EQ(std::complex<double>(rank_avg, d).real(), f(1).real());

  EXPECT_DOUBLE_EQ(err_expected(0).real(), err_overwriting(0).real());
  EXPECT_DOUBLE_EQ(err_expected(0).imag(), err_overwriting(0).imag());
  EXPECT_DOUBLE_EQ(err_expected(1).real(), err_overwriting(1).real());
  EXPECT_DOUBLE_EQ(err_expected(1).imag(), err_overwriting(1).imag());
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceScalar) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<FunctionDomain, FunctionDomain>;

  dca::func::function<double, FunctionDomain> f("f");
  dca::func::function<double, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); ++i) {
    f(i) = i * rank_;
    f_mean(i) = f(i);
  }
  sum_interface_.sum(f_mean);
  f_mean /= size_;

  dca::func::function<double, CovarianceDomain> covariance_expected("expected");
  covariance_expected(0, 0) = 0.;
  covariance_expected(0, 1) = 0.;
  covariance_expected(0, 2) = 0.;
  covariance_expected(0, 3) = 0.;
  covariance_expected(1, 1) = 5.25;
  covariance_expected(1, 2) = 10.5;
  covariance_expected(1, 3) = 15.75;
  covariance_expected(2, 2) = 21.0;
  covariance_expected(2, 3) = 31.5;
  covariance_expected(3, 3) = 47.25;

  for (int i = 0; i < f.size(); ++i)
    for (int j = i + 1; j < f.size(); ++j)
      covariance_expected(j, i) = covariance_expected(i, j);

  dca::func::function<double, CovarianceDomain> covariance_from_computeCovariance(
      "computeCovariance");
  dca::func::function<double, CovarianceDomain> covariance_from_computeCovarianceAndAvg(
      "computeCovarianceAndAvg");

  // Compute covariance matrix with respect to the *precomputed* mean of f.
  sum_interface_.computeCovariance(f, f_mean, covariance_from_computeCovariance);
  // Compute covariance matrix with respect to the mean of f.
  sum_interface_.computeCovarianceAndAvg(f, covariance_from_computeCovarianceAndAvg);

  for (int i = 0; i < f_mean.size(); ++i)
    EXPECT_DOUBLE_EQ(f_mean(i), f(i));

  // The covariance matrices are identical as both are computed with respect to the mean of f.
  for (int i = 0; i < covariance_expected.size(); ++i) {
    EXPECT_DOUBLE_EQ(covariance_expected(i), covariance_from_computeCovariance(i));
    EXPECT_DOUBLE_EQ(covariance_expected(i), covariance_from_computeCovarianceAndAvg(i));
  }
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceComplex) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<dca::func::dmn_0<dca::func::dmn<4, int>>,
                                                   dca::func::dmn_0<dca::func::dmn<4, int>>>;

  dca::func::function<std::complex<double>, FunctionDomain> f("f");
  dca::func::function<std::complex<double>, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); ++i) {
    f(i) = std::complex<double>(rank_ * i, rank_ * (i + f.size()));
    f_mean(i) = f(i);
  }
  sum_interface_.sum(f_mean);
  f_mean /= size_;

  dca::func::function<double, CovarianceDomain> covariance_expected("expected");
  covariance_expected(0, 0) = 0.;
  covariance_expected(0, 1) = 0.;
  covariance_expected(0, 2) = 0.;
  covariance_expected(0, 3) = 0.;
  covariance_expected(1, 1) = 5.25;
  covariance_expected(1, 2) = 10.5;
  covariance_expected(1, 3) = 15.75;
  covariance_expected(2, 2) = 21.0;
  covariance_expected(2, 3) = 31.5;
  covariance_expected(3, 3) = 47.25;

  for (int i = 0; i < 2 * f.size(); ++i)
    for (int j = i + 1; j < 2 * f.size(); ++j)
      covariance_expected(j, i) = covariance_expected(i, j);

  dca::func::function<double, CovarianceDomain> covariance_from_computeCovariance(
      "computeCovariance");
  dca::func::function<double, CovarianceDomain> covariance_from_computeCovarianceAndAvg(
      "computeCovarianceAndAvg");

  // Compute covariance matrix with respect to the *precomputed* mean of f.
  sum_interface_.computeCovariance(f, f_mean, covariance_from_computeCovariance);
  // Compute covariance matrix with respect to the mean of f.
  sum_interface_.computeCovarianceAndAvg(f, covariance_from_computeCovarianceAndAvg);

  for (int i = 0; i < f_mean.size(); ++i) {
    EXPECT_DOUBLE_EQ(f_mean(i).real(), f(i).real());
    EXPECT_DOUBLE_EQ(f_mean(i).imag(), f(i).imag());
  }

  // The covariance matrices are identical as both are computed with respect to the mean of f.
  for (int i = 0; i < covariance_expected.size(); ++i) {
    EXPECT_DOUBLE_EQ(covariance_expected(i), covariance_from_computeCovariance(i));
    EXPECT_DOUBLE_EQ(covariance_expected(i), covariance_from_computeCovarianceAndAvg(i));
  }
}

TEST_F(MPICollectiveSumTest, AvgNormalizedMomenta) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  const std::vector<int> orders{3, 4};

  dca::func::function<double, FunctionDomain> f("f");
  std::vector<double> momenta(orders.size());

  for (int i = 0; i < f.size(); i++) {
    f(i) = rank_ * (i + 1);
  }

  momenta = sum_interface_.avgNormalizedMomenta(f, orders);
  // Expected values obtained with python
  EXPECT_NEAR(0., momenta[0], 1e-8);
  EXPECT_NEAR(1.76190476, momenta[1], 1e-8);
}

TEST_F(MPICollectiveSumTest, DelayedSum) {
  using TestDomain1 = dca::func::dmn_0<dca::func::dmn<2, char>>;
  using TestDomain2 = dca::func::dmn_0<dca::func::dmn<3, char>>;

  dca::func::function<double, TestDomain1> function_test("test");
  dca::func::function<double, TestDomain1> function_expected("expected");

  dca::func::function<std::complex<double>, TestDomain2> function_cmplx_test("test");
  dca::func::function<std::complex<double>, TestDomain2> function_cmplx_expected("expected");

  for (int i = 0; i < function_test.size(); i++) {
    function_test(i) = i * rank_;
    function_cmplx_test(i) = std::complex<double>(i * rank_, 1);
    function_cmplx_expected(i) = std::complex<double>(i * size_ * (size_ - 1) / 2, size_);
    function_expected(i) = i * size_ * (size_ - 1) / 2;
  }

  double scalar_test = rank_;
  double scalar_expected = size_ * (size_ - 1) / 2;

  sum_interface_.delayedSum(function_test);
  sum_interface_.delayedSum(scalar_test);
  sum_interface_.delayedSum(function_cmplx_test);

  sum_interface_.resolveSums();

  EXPECT_EQ(scalar_expected, scalar_test);
  EXPECT_EQ(function_expected, function_test);
  EXPECT_EQ(function_cmplx_expected, function_cmplx_test);
}

int main(int argc, char** argv) {
  int result = 0;

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
