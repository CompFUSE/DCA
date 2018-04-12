// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests richardson_lucy_deconvolution.hpp

#include "dca/math/inference/richardson_lucy_deconvolution.hpp"

#include <utility>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"

class RichardsonLucyDeconvolutionTest : public ::testing::Test {
protected:
  using ClusterDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using HostDmn = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using OtherDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

  RichardsonLucyDeconvolutionTest()
      : p_cluster_(std::make_pair(ClusterDmn::dmn_size(), HostDmn::dmn_size()),
                   "projection-operator-cluster"),
        p_host_(HostDmn::dmn_size(), "projection-operator-host") {}

  dca::linalg::Matrix<double, dca::linalg::CPU> p_cluster_;
  dca::linalg::Matrix<double, dca::linalg::CPU> p_host_;

  dca::func::function<double, dca::func::dmn_variadic<ClusterDmn, OtherDmn>> source_;
  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> source_interpolated_;
  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> target_;
  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> target_convoluted_;
};

TEST_F(RichardsonLucyDeconvolutionTest, IdentityProjectionOperator) {
  // Projection operator = identity matrix.
  p_cluster_(0, 0) = p_cluster_(1, 1) = 1.;
  p_host_(0, 0) = p_host_(1, 1) = p_host_(2, 2) = p_host_(3, 3) = 1.;

  const double tolerance = 1.e-3;
  const int max_iterations = 3;

  dca::math::inference::RichardsonLucyDeconvolution<ClusterDmn, HostDmn, OtherDmn> deconvolution(
      p_cluster_, p_host_, tolerance, max_iterations);

  // Trivial function not close to or crossing zero.
  // Do not use 1., since this would be equal to the initial guess.
  source_ = 2.;
  source_interpolated_ = 2.;

  //
  // Test findTargetFunction(source, source_interpolated, target).
  //
  int iterations = deconvolution.findTargetFunction(source_, source_interpolated_, target_);

  EXPECT_EQ(1, iterations);

  for (int i = 0; i < source_interpolated_.size(); ++i) {
    EXPECT_DOUBLE_EQ(source_interpolated_(i), target_(i));
  }

  //
  // Test findTargetFunction(source, source_interpolated, target, target_convoluted).
  //
  target_ = 0.;

  iterations =
      deconvolution.findTargetFunction(source_, source_interpolated_, target_, target_convoluted_);

  EXPECT_EQ(1, iterations);

  for (int i = 0; i < source_interpolated_.size(); ++i) {
    EXPECT_DOUBLE_EQ(source_interpolated_(i), target_(i));
    EXPECT_DOUBLE_EQ(source_interpolated_(i), target_convoluted_(i));
  }
}

TEST_F(RichardsonLucyDeconvolutionTest, SourceCloseToZero) {
  const double tolerance = 1.e-3;
  const int max_iterations = 3;

  dca::math::inference::RichardsonLucyDeconvolution<ClusterDmn, HostDmn, OtherDmn> deconvolution(
      p_cluster_, p_host_, tolerance, max_iterations);

  source_interpolated_ = 1.e-8;

  EXPECT_THROW(deconvolution.findTargetFunction(source_, source_interpolated_, target_),
               std::invalid_argument);
}
