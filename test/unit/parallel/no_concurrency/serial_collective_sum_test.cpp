// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests serial_collective_sum.hpp.

#include "dca/parallel/no_concurrency/serial_collective_sum.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

class SerialCollectiveSumTest : public ::testing::Test {
protected:
  using Domain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  using RealFunction = dca::func::function<double, Domain>;
  using ComplexFunction = dca::func::function<std::complex<double>, Domain>;

  SerialCollectiveSumTest() : d0_(3.14), d1_(2.72), d2_(42.) {}

  virtual void SetUp() {
    f_real_(0) = d0_;
    f_real_(1) = d1_;
    f_real_(2) = d2_;

    f_complex_(0) = std::complex<double>(d0_, d1_);
    f_complex_(1) = std::complex<double>(d1_, d2_);
    f_complex_(2) = std::complex<double>(d2_, d0_);
  }

  dca::parallel::SerialCollectiveSum sum_interface_;

  RealFunction f_real_;
  ComplexFunction f_complex_;

  const double d0_;
  const double d1_;
  const double d2_;
};

TEST_F(SerialCollectiveSumTest, Sum) {
  // Scalar
  double s = d0_;
  sum_interface_.sum(s);
  EXPECT_EQ(d0_, s);

  // std::vector
  const std::vector<int> v_check{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<int> v(v_check);
  sum_interface_.sum(v);
  EXPECT_EQ(v_check, v);

  // std::map<std::string, std::vector>

  // dca::func::function<Scalar, Domain>
  sum_interface_.sum(f_real_);

  EXPECT_EQ(d0_, f_real_(0));
  EXPECT_EQ(d1_, f_real_(1));
  EXPECT_EQ(d2_, f_real_(2));

  // dca::func::function<Scalar, Domain> (separate f_real_in and f_real_out)

  // dca::func::function<std::vector, Domain>

  // dca::linalg::Vector

  // dca::linalg::Matrix
}

TEST_F(SerialCollectiveSumTest, SumAndAverage) {
  const int measurements = 2000;

  // Scalar
  double s = d0_;
  sum_interface_.sum_and_average(s, measurements);
  EXPECT_EQ(d0_ / measurements, s);

  // dca::func::function
  sum_interface_.sum_and_average(f_real_, measurements);

  EXPECT_EQ(d0_ / measurements, f_real_(0));
  EXPECT_EQ(d1_ / measurements, f_real_(1));
  EXPECT_EQ(d2_ / measurements, f_real_(2));
}

TEST_F(SerialCollectiveSumTest, AverageAndComputeStddev) {
  RealFunction f_real_stddev;

  sum_interface_.average_and_compute_stddev(f_real_, f_real_stddev);

  EXPECT_EQ(d0_, f_real_(0));
  EXPECT_EQ(d1_, f_real_(1));
  EXPECT_EQ(d2_, f_real_(2));

  EXPECT_EQ(0., f_real_stddev(0));
  EXPECT_EQ(0., f_real_stddev(1));
  EXPECT_EQ(0., f_real_stddev(2));
}

TEST_F(SerialCollectiveSumTest, LeaveOneOutAvg) {
  // Scalar
  double s = d0_;
  sum_interface_.leaveOneOutAvg(s);
  EXPECT_EQ(d0_, s);

  // dca::func::function
  sum_interface_.leaveOneOutAvg(f_real_);
  EXPECT_EQ(d0_, f_real_(0));
  EXPECT_EQ(d1_, f_real_(1));
  EXPECT_EQ(d2_, f_real_(2));
}

TEST_F(SerialCollectiveSumTest, JackknifeErrorReal) {
  const RealFunction f_real_copy(f_real_);

  // Pass true.
  auto err_true = sum_interface_.jackknifeError(f_real_, true);

  for (int i = 0; i < f_real_.size(); ++i) {
    EXPECT_EQ(f_real_copy(i), f_real_(i));
    EXPECT_EQ(0., err_true(i));
  }

  // Pass false.
  auto err_false = sum_interface_.jackknifeError(f_real_, false);

  for (int i = 0; i < f_real_.size(); ++i) {
    EXPECT_EQ(f_real_copy(i), f_real_(i));
    EXPECT_EQ(0., err_false(i));
  }
}

TEST_F(SerialCollectiveSumTest, JackknifeErrorComplex) {
  const ComplexFunction f_complex_copy(f_complex_);

  // Pass true.
  auto err_true = sum_interface_.jackknifeError(f_complex_, true);

  for (int i = 0; i < f_complex_.size(); ++i) {
    EXPECT_EQ(f_complex_copy(i), f_complex_(i));
    EXPECT_EQ(0., err_true(i));
  }

  // Pass false.
  auto err_false = sum_interface_.jackknifeError(f_complex_, false);

  for (int i = 0; i < f_complex_.size(); ++i) {
    EXPECT_EQ(f_complex_copy(i), f_complex_(i));
    EXPECT_EQ(0., err_false(i));
  }
}
