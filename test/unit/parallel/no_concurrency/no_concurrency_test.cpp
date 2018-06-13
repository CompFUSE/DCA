// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests no_concurrency.hpp.

#include "dca/parallel/no_concurrency/no_concurrency.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

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
  using Domain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  Domain dmn;

  std::pair<int, int> bounds = concurrency_.get_bounds(dmn);
  EXPECT_EQ(0, bounds.first);
  EXPECT_EQ(3, bounds.second);
}
