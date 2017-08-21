// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests some general behaviour of the random number generators.

#include "dca/math/random/random.hpp"
#include <utility>  // for std::move
#include <vector>
#include "gtest/gtest.h"

template <typename Generator>
class RandomTest : public ::testing::Test {};

using Generators = ::testing::Types<dca::math::random::StdRandomWrapper<std::mt19937_64>>;

TYPED_TEST_CASE(RandomTest, Generators);

TYPED_TEST(RandomTest, MoveConstruction) {
  using RngType = TypeParam;

  RngType rng_1(0, 1);
  const int global_id = rng_1.getGlobalId();
  const uint64_t seed = rng_1.getSeed();

  RngType rng_2(std::move(rng_1));
  EXPECT_EQ(global_id, rng_2.getGlobalId());
  EXPECT_EQ(seed, rng_2.getSeed());
}
