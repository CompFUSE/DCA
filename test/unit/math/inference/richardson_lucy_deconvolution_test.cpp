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

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"

TEST(RichardsonLucyDeconvolutionTest, IdentityProjectionOperator) {
  using DeconvolutionDmn = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using OtherDmn = dca::func::dmn_0<dca::func::dmn<1, int>>;

  const double tolerance = 1.e-3;
  const int max_iterations = 3;

  dca::math::inference::RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn> deconvolution(
      tolerance, max_iterations);

  // Projection operator = identity matrix.
  dca::linalg::Matrix<double, dca::linalg::CPU> p(4, "projection-operator");
  p(0, 0) = p(1, 1) = p(2, 2) = p(3, 3) = 1.;

  // Trivial function not crossing zero.
  dca::func::function<double, dca::func::dmn_variadic<DeconvolutionDmn, OtherDmn>> source;
  for (int i = 0; i < source.size(); ++i) {
    source(i) = 1.;
  }

  dca::func::function<double, dca::func::dmn_variadic<DeconvolutionDmn, OtherDmn>> target;

  const int iterations = deconvolution.execute(p, source, target);

  EXPECT_EQ(1, iterations);

  for (int i = 0; i < source.size(); ++i) {
    EXPECT_DOUBLE_EQ(source(i), target(i));
  }
}
