// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests ignore_unused.hpp

#include "dca/util/ignore_unused.hpp"
#include "gtest/gtest.h"

TEST(IgnoreUnusedTest, basicTest) {
  int i = 42;
  double d = 3.14;

  dca::util::ignoreUnused();
  dca::util::ignoreUnused(i);
  dca::util::ignoreUnused(d);
}
