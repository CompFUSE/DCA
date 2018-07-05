// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the ceilDiv function.

#include "dca/util/integer_division.hpp"
#include <vector>
#include "gtest/gtest.h"

template <typename ScalarType>
class IntDivTest : public ::testing::Test {};

typedef ::testing::Types<short, int, long, long long, unsigned short, unsigned int, unsigned long,
                         unsigned long long>
    IntegerTypes;
TYPED_TEST_CASE(IntDivTest, IntegerTypes);

TYPED_TEST(IntDivTest, CeilDiv) {
  using IntType = TypeParam;
  std::vector<IntType> a = {0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  std::vector<IntType> res = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3};

  IntType b(5);
  for (int i = 0; i < a.size(); ++i)
    EXPECT_EQ(res[i], dca::util::ceilDiv(a[i], b));
}
