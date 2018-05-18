// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the file pack_operations.hpp

#include "dca/util/pack_operations.hpp"

#include "gtest/gtest.h"

TEST(PackOperationsTest, IfAll) {
  constexpr bool b1 = dca::util::if_all<true, std::is_integral<int>::value, 1>::value;
  EXPECT_TRUE(b1);

  constexpr bool b2 = dca::util::if_all<true, std::is_integral<double>::value, 1>::value;
  EXPECT_FALSE(b2);
}

TEST(PackOperationsTest, Product) {
  const int result = dca::util::product(3, 2, 4);
  EXPECT_EQ(3 * 2 * 4, result);
}
