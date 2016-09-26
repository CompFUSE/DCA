// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests no_concurrency.hpp.

#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "gtest/gtest.h"
#include "comp_library/function_library/domains/special_domains/dmn.h"
#include "comp_library/function_library/domains/special_domains/dmn_0.h"

class NoConcurrencyTest : public ::testing::Test {
protected:
  NoConcurrencyTest() : concurrency_(0, nullptr) {}
  dca::parallel::NoConcurrency concurrency_;
};

TEST_F(NoConcurrencyTest, IDsAndSize) {
  EXPECT_EQ(0, concurrency_.id());
  EXPECT_EQ(1, concurrency_.number_of_processors());
  EXPECT_EQ(0, concurrency_.first());
  EXPECT_EQ(0, concurrency_.last());
}

TEST_F(NoConcurrencyTest, Broadcast) {
  const int root_id = 0;

  double d = 3.14;
  EXPECT_TRUE(concurrency_.broadcast(d, root_id));

  unsigned int u = 123;
  EXPECT_TRUE(concurrency_.broadcast(u));

  int i = -42;
  EXPECT_TRUE(concurrency_.broadcast_object(i, root_id));

  float f = 2.72;
  EXPECT_TRUE(concurrency_.broadcast_object(f));
}

TEST_F(NoConcurrencyTest, GetBounds) {
  using Domain = dmn_0<dmn<3, int>>;
  Domain dmn;

  std::pair<int, int> bounds = concurrency_.get_bounds(dmn);
  EXPECT_EQ(0, bounds.first);
  EXPECT_EQ(3, bounds.second);
}
