// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the random number library utility functions.

#include "dca/math/random/random_utils.hpp"
#include <algorithm>  // for std::sort and std::unique
#include <cstdint>    // for uint64_t
#include <numeric>    // for std::iota
#include <vector>
#include "gtest/gtest.h"

using dca::math::random::detail::getGlobalId;
using dca::math::random::detail::generateSeed;
using dca::math::random::detail::hash;

TEST(GetGlobalIdTest, FirstId) {
  EXPECT_EQ(0, getGlobalId(0, 0, 1));
  EXPECT_EQ(0, getGlobalId(0, 0, 10));
}

TEST(GetGlobalIdTest, IncreaseProcId) {
  EXPECT_EQ(1, getGlobalId(0, 1, 4));
  EXPECT_EQ(2, getGlobalId(0, 2, 25));
}

TEST(GetGlobalIdTest, IncreaseLocalId) {
  EXPECT_EQ(4, getGlobalId(1, 0, 4));
  EXPECT_EQ(50, getGlobalId(2, 0, 25));
}

TEST(GetGlobalIdTest, Sequence) {
  int num_procs = 10;
  int num_threads = 4;

  std::vector<int> generated_ids;
  std::vector<int> expected_ids(num_procs * num_threads);
  std::iota(expected_ids.begin(), expected_ids.end(), 0);

  for (int local_id = 0; local_id < num_threads; ++local_id) {
    for (int proc_id = 0; proc_id < num_procs; ++proc_id) {
      generated_ids.push_back(getGlobalId(local_id, proc_id, num_procs));
    }
  }

  EXPECT_EQ(expected_ids, generated_ids);
}

TEST(generateSeedTest, DefaultSeed) {
  EXPECT_EQ(0, generateSeed(0));
  EXPECT_EQ(17886235652272626782u, generateSeed(100));

  EXPECT_EQ(generateSeed(100), generateSeed(100, 0));
}

TEST(generateSeedTest, CustomSeed) {
  EXPECT_EQ(3924199218512185873, generateSeed(0, 42));
  EXPECT_EQ(13794760588299155443u, generateSeed(1234, 42));
  EXPECT_EQ(5039178850028241505, generateSeed(25, -137));

  EXPECT_EQ(generateSeed(142), generateSeed(100, 42));

  EXPECT_NE(generateSeed(1234), generateSeed(1234, 42));
  EXPECT_NE(generateSeed(1234, 42), generateSeed(1234, 137));
}

TEST(generateSeedTest, Uniqueness) {
  // Check for seed collisions when generating 10^7 seeds.
  int max_global_id = 10000000;
  std::vector<uint64_t> seeds;

  for (int global_id = 0; global_id < max_global_id; ++global_id)
    seeds.push_back(generateSeed(global_id, 0));

  std::sort(seeds.begin(), seeds.end());
  std::vector<uint64_t>::iterator new_end = std::unique(seeds.begin(), seeds.end());

  EXPECT_EQ(max_global_id, new_end - seeds.begin());
}

TEST(HashTest, Zero) {
  EXPECT_EQ(0, hash(0));
}

TEST(HashTest, Positive) {
  EXPECT_EQ(3924199218512185873, hash(42));
  EXPECT_EQ(4496399509612857899, hash(22737172172));
}

TEST(HashTest, Negative) {
  int i = -137;

  EXPECT_EQ(18446744073709551479u, uint64_t(i));
  EXPECT_EQ(hash(i), hash(18446744073709551479u));
  EXPECT_EQ(13902458719468866193u, hash(i));
}
