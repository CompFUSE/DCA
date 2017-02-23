// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests comparison_methods.hpp.

#include "dca/math/util/comparison_methods.hpp"
#include <complex>
#include <utility>
#include "gtest/gtest.h"

TEST(ComparisonMethodsTest, SusceptibilityPairLess) {
  // Real
  const std::pair<double, int> x1(0.9, 42);
  const std::pair<double, int> y1(1.2, 41);
  EXPECT_TRUE(dca::math::util::susceptibilityPairLess(x1, y1));
  EXPECT_FALSE(dca::math::util::susceptibilityPairLess(y1, x1));

  const std::pair<int, double> x2(0, 4.2);
  const std::pair<int, double> y2(3, 4.1);
  EXPECT_TRUE(dca::math::util::susceptibilityPairLess(x2, y2));

  const std::pair<double, int> x3(1., 99);
  const std::pair<double, int> y3(1., 100);
  EXPECT_FALSE(dca::math::util::susceptibilityPairLess(x3, y3));
}
