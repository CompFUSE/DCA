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

TEST(RichardsonLucyDeconvolutionTest, IdentityProjectionOperator) {
  using ClusterDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using HostDmn = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using OtherDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

  // Projection operator = identity matrix.
  dca::linalg::Matrix<double, dca::linalg::CPU> p_cluster(
      std::make_pair(ClusterDmn::dmn_size(), HostDmn::dmn_size()), "projection-operator-cluster");
  p_cluster(0, 0) = p_cluster(1, 1) = 1.;

  dca::linalg::Matrix<double, dca::linalg::CPU> p_host(HostDmn::dmn_size(),
                                                       "projection-operator-host");
  p_host(0, 0) = p_host(1, 1) = p_host(2, 2) = p_host(3, 3) = 1.;

  const double tolerance = 1.e-3;
  const int max_iterations = 3;

  dca::math::inference::RichardsonLucyDeconvolution<ClusterDmn, HostDmn, OtherDmn> deconvolution(
      p_cluster, p_host, tolerance, max_iterations);

  // Trivial function not crossing zero.
  dca::func::function<double, dca::func::dmn_variadic<ClusterDmn, OtherDmn>> source;
  // Do not use 1., since this would be equal to the initial guess.
  source = 2.;

  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> source_interpolated;
  source_interpolated = 2.;

  //
  // Test findTargetFunction(source, source_interpolated, target).
  //
  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> target;

  int iterations = deconvolution.findTargetFunction(source, source_interpolated, target);

  EXPECT_EQ(1, iterations);

  for (int i = 0; i < source_interpolated.size(); ++i) {
    EXPECT_DOUBLE_EQ(source_interpolated(i), target(i));
  }

  //
  // Test findTargetFunction(source, source_interpolated, target, target_convoluted).
  //
  target = 0.;
  dca::func::function<double, dca::func::dmn_variadic<HostDmn, OtherDmn>> target_convoluted;

  iterations =
      deconvolution.findTargetFunction(source, source_interpolated, target, target_convoluted);

  EXPECT_EQ(1, iterations);

  for (int i = 0; i < source_interpolated.size(); ++i) {
    EXPECT_DOUBLE_EQ(source_interpolated(i), target(i));
    EXPECT_DOUBLE_EQ(source_interpolated(i), target_convoluted(i));
  }
}
