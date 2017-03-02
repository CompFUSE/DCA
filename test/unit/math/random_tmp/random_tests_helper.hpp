// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides helper functions for the random number library unit tests.

#ifndef DCA_TEST_UNIT_MATH_RANDOM_RANDOM_TESTS_HELPER_HPP
#define DCA_TEST_UNIT_MATH_RANDOM_RANDOM_TESTS_HELPER_HPP

#include <cassert>
#include "gtest/gtest.h"

namespace dca {
namespace testing {
// dca::testing::

// Checks whether the drawn pseudo-random numbers are in the interval [0, 1).
template <typename RandomNumberGenerator>
void inUnitInterval(RandomNumberGenerator& rng, const int draws) {
  assert(draws > 0);

  for (int i = 0; i < draws; ++i) {
    double rn = rng();
    EXPECT_TRUE(rn >= 0.);
    EXPECT_TRUE(rn < 1.);
  }
}

}  // testing
}  // dca

#endif  // DCA_TEST_UNIT_MATH_RANDOM_RANDOM_TESTS_HELPER_HPP
