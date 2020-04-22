// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the file pack_operations.hpp

#include "dca/util/pack_operations.hpp"

#include "gtest/gtest.h"

TEST(PackOperationsTest, IfAll) {
  constexpr bool b1 = dca::util::ifAll(true, std::is_integral_v<int>, 1);
  EXPECT_TRUE(b1);

  constexpr bool b2 = dca::util::ifAll(true, std::is_integral_v<double>, 1);
  EXPECT_FALSE(b2);
}

TEST(PackOperationsTest, Product) {
  const int result = dca::util::product(3, 2, 4);
  EXPECT_EQ(3 * 2 * 4, result);
}

TEST(PackOperationsTest, SizeSum) {
  using dca::util::size_sum;
  EXPECT_EQ(sizeof(char), size_sum<char>);

  constexpr unsigned sum = size_sum<int, double, unsigned>;
  EXPECT_EQ(sizeof(int) + sizeof(double) + sizeof(unsigned), sum);

  using List = dca::util::Typelist<int, double, unsigned>;
  EXPECT_EQ(sum, size_sum<List>);

  using EmptyList = dca::util::Typelist<>;
  EXPECT_EQ(0, size_sum<EmptyList>);
}
