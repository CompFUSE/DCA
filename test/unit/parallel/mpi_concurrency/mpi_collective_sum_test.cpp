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

TEST_F(MPICollectiveSumTest, LeaveOneOutAverage) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  dca::func::function<double, TestDomain> f;

  std::vector<double>values(size_);
  double sum = 0.;
  for(int i =0 ;i < size_; i++) {
    values[i] = 10. + i;
    sum += values[i];
  }
  f = values[rank_];
  double scalar = values[rank_];

  sum_interface_.leaveOneOutAvg(f);
  sum_interface_.leaveOneOutAvg(scalar);
  double expected = (sum - values[rank_]) / (size_-1);
  for(int i =0; i < f.size(); i++)
    EXPECT_NEAR(expected, f(i), 1e-10);
  EXPECT_NEAR(expected, scalar, 1e-10);
}

TEST_F(MPICollectiveSumTest, JackKnifeTestReal) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  dca::func::function<double, TestDomain> f;

  auto val = [](int rank){return rank + 10.;};
  f = val(rank_);

  sum_interface_.leaveOneOutAvg(f);

  for(int i =0; i < f.size(); i++)
    f(i) = std::pow(f(i), 3);

  auto err = sum_interface_.jackKnifeError(f);

  // Expected error computed with python.
  const double expected = 473.7688;
  for(int i =0; i < err.size(); i++)
    EXPECT_NEAR(expected, err(i), 1e-3);


}

TEST_F(MPICollectiveSumTest, JackKnifeTestComplex){
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  dca::func::function<std::complex<double>, TestDomain> f_c;
  auto val = [](int rank){return rank + 10.;};

  for(int i =0; i < f_c.size(); i++) {
    f_c(i).real(val(rank_));
    f_c(i).imag(val(rank_));
  }
  sum_interface_.leaveOneOutAvg(f_c);
  for(int i =0; i < f_c.size(); i++) {
    f_c(i).real(std::pow(std::real(f_c(i)), 3));
    f_c(i).imag(std::pow(std::imag(f_c(i)), 3));
  }
  auto err_c = sum_interface_.jackKnifeError(f_c);

  const double expected = 473.7688;
  for(int i =0; i < err_c.size(); i++) {
    EXPECT_NEAR(expected, std::real(err_c(i)), 1e-3);
    EXPECT_NEAR(expected, std::imag(err_c(i)), 1e-3);
  }
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceScalar) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<FunctionDomain, FunctionDomain>;

  dca::func::function<double, CovarianceDomain> covariance_test("test");
  dca::func::function<double, CovarianceDomain> covariance_expected("expected");

  dca::func::function<double, FunctionDomain> f("f");
  dca::func::function<double, FunctionDomain> f_mean("f_mean");

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
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<dca::func::dmn_0<dca::func::dmn<4, int>>,
                                                   dca::func::dmn_0<dca::func::dmn<4, int>>>;

  dca::func::function<double, CovarianceDomain> covariance_test("test");
  dca::func::function<double, CovarianceDomain> covariance_expected("expected");

  dca::func::function<std::complex<double>, FunctionDomain> f("f");
  dca::func::function<std::complex<double>, FunctionDomain> f_mean("f_mean");

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
