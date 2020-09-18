// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the new methods in type_lists.hpp

#include "dca/util/type_list.hpp"

#include "gtest/gtest.h"

TEST(TypeListTest, Sublist) {
  using List = dca::util::Typelist<int, float, double, char>;
  using Sublist = dca::util::Sublist<2, List>;

  EXPECT_EQ(4, dca::util::Length<List>::value);
  EXPECT_EQ(2, dca::util::Length<Sublist>::value);

  constexpr bool elements_eq = std::is_same<dca::util::Typelist<int, float>, Sublist>::value;
  EXPECT_TRUE(elements_eq);
}
