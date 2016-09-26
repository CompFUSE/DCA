// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Andrei Plamada (plamada@phys.ethz.ch)
//
// This file tests mpi_collective_sum.hpp.

#include "dca/parallel/mpi_concurrency/mpi_collective_sum.hpp"
#include <complex>
#include "gtest/gtest.h"
#include "dca/testing/minimalist_printer.hpp"

class MPICollectiveSumTest : public ::testing::Test {
protected:
  MPICollectiveSumTest() : sum_interface_(grouping_) {
    grouping_.set();
    size_ = grouping_.get_Nr_threads();
    rank_ = grouping_.get_id();
  }

  int size_;
  int rank_;

  dca::parallel::MPIProcessorGrouping grouping_;
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
  using TestDomain = dmn_0<dmn<2, int>>;

  FUNC_LIB::function<double, TestDomain> function_test("test");
  FUNC_LIB::function<double, TestDomain> function_expected("expected");

  for (int i = 0; i < function_test.size(); i++)
    function_test(i) = i * rank_;

  sum_interface_.sum(function_test);

  for (int i = 0; i < function_test.size(); i++)
    function_expected(i) = i * size_ * (size_ - 1) / 2;

  for (int i = 0; i < function_test.size(); i++)
    EXPECT_EQ(function_expected(i), function_test(i));
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceScalar) {
  using FunctionDomain = dmn_0<dmn<4, int>>;
  using CovarianceDomain = dmn_variadic<FunctionDomain, FunctionDomain>;

  FUNC_LIB::function<double, CovarianceDomain> covariance_test("test");
  FUNC_LIB::function<double, CovarianceDomain> covariance_expected("expected");

  FUNC_LIB::function<double, FunctionDomain> f("f");
  FUNC_LIB::function<double, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); i++) {
    f(i) = i * rank_;
    f_mean(i) = f(i);
  }
  sum_interface_.sum(f_mean);
  f_mean /= size_;  // f_mean contains the mean of f

  sum_interface_.computeCovariance(f, f_mean, covariance_test);

  covariance_expected(0, 0) = covariance_expected(0, 1) = covariance_expected(0, 2) =
      covariance_expected(0, 3) = 0;
  covariance_expected(1, 1) = 5.25;
  covariance_expected(1, 2) = 10.5;
  covariance_expected(1, 3) = 15.75;
  covariance_expected(2, 2) = 21.0;
  covariance_expected(2, 3) = 31.5;
  covariance_expected(3, 3) = 47.25;

  for (int i = 0; i < f.size(); i++)
    for (int j = i + 1; j < f.size(); j++)
      covariance_expected(j, i) = covariance_expected(i, j);

  for (int i = 0; i < covariance_test.size(); i++)
    EXPECT_EQ(covariance_expected(i), covariance_test(i));
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceComplex) {
  using FunctionDomain = dmn_0<dmn<2, int>>;
  using CovarianceDomain = dmn_variadic<dmn_0<dmn<4, int>>, dmn_0<dmn<4, int>>>;

  FUNC_LIB::function<double, CovarianceDomain> covariance_test("test");
  FUNC_LIB::function<double, CovarianceDomain> covariance_expected("expected");

  FUNC_LIB::function<std::complex<double>, FunctionDomain> f("f");
  FUNC_LIB::function<std::complex<double>, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); i++) {
    f(i) = std::complex<double>(rank_ * i, rank_ * (i + f.size()));
    f_mean(i) = f(i);
  }

  sum_interface_.sum(f_mean);
  f_mean /= size_;  // f_mean contains the mean of f

  sum_interface_.computeCovariance(f, f_mean, covariance_test);

  covariance_expected(0, 0) = covariance_expected(0, 1) = covariance_expected(0, 2) =
      covariance_expected(0, 3) = 0;
  covariance_expected(1, 1) = 5.25;
  covariance_expected(1, 2) = 10.5;
  covariance_expected(1, 3) = 15.75;
  covariance_expected(2, 2) = 21.0;
  covariance_expected(2, 3) = 31.5;
  covariance_expected(3, 3) = 47.25;

  for (int i = 0; i < 2 * f.size(); i++)
    for (int j = i + 1; j < 2 * f.size(); j++)
      covariance_expected(j, i) = covariance_expected(i, j);

  for (int i = 0; i < covariance_test.size(); i++)
    EXPECT_EQ(covariance_expected(i), covariance_test(i));
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
