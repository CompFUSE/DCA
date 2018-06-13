// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file tests difference.hpp.

#include "dca/function/util/difference.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

template <int num_elements, typename Scalar = double>
using Function = dca::func::function<Scalar, dca::func::dmn_0<dca::func::dmn<num_elements>>>;

TEST(DifferenceTest, Double) {
  Function<4> f1, f2;
  f1 = 1.;

  f2(0) = 1.3;
  f2(1) = 0.7;
  f2(2) = 2;
  f2(3) = 1.4;

  dca::func::util::Difference diff = dca::func::util::difference(f1, f2);

  EXPECT_NEAR(0.5, diff.l1, 1e-4);
  EXPECT_NEAR(0.57879, diff.l2, 1e-4);
  EXPECT_NEAR(1, diff.l_inf, 1e-4);
}
