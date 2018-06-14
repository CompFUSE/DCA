// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the construction of random number generators in a std::vector.
// The test is put in a separate file (compilation unit) to keep control over static variables.

#include "dca/math/random/random.hpp"
#include <vector>
#include "gtest/gtest.h"

template <typename Generator>
class RandomTest : public ::testing::Test {};

using Generators = ::testing::Types<dca::math::random::StdRandomWrapper<std::mt19937_64>>;

TYPED_TEST_CASE(RandomTest, Generators);

TYPED_TEST(RandomTest, RngVectorConstruction) {
  using RngType = TypeParam;

  const int num_rngs = 5;
  const int proc_id = 1;
  const int num_procs = 4;
  const int seed = 42;

  std::vector<RngType> rngs;

  for (int i = 0; i < num_rngs; ++i)
    rngs.emplace_back(proc_id, num_procs, seed);

  std::vector<int> expected_ids{1, 5, 9, 13, 17};

  for (int i = 0; i < rngs.size(); ++i)
    EXPECT_EQ(expected_ids[i], rngs[i].getGlobalId());
}
