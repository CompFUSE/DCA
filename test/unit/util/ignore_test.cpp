// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests ignore.hpp

#include <vector>
#include "dca/util/ignore.hpp"
#include "gtest/gtest.h"

TEST(IgnoreUnusedTest, NoArguments) {
  dca::util::ignoreUnused();
}

TEST(IgnoreUnusedTest, FundamentalTypes) {
  int i = 42;
  double d = 3.14;
  dca::util::ignoreUnused(i);
  dca::util::ignoreUnused(d);
}

TEST(IgnoreUnusedTest, StdTypes) {
  std::vector<double> v(10, 3.14);
  dca::util::ignoreUnused(v);
}

class MyClass {};

TEST(IgnoreUnusedTest, CustomTypes) {
  MyClass c;
  dca::util::ignoreUnused(c);
}

TEST(IgnoreUnusedTest, MultipleArguments) {
  float f = 1.2f;
  bool b = true;
  dca::util::ignoreUnused(f, b);
}
