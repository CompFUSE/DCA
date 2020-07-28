// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the RandomAccessSet class.

#include "dca/util/containers/random_access_set.hpp"

#include <set>
#include <random>
#include <string>

#include "gtest/gtest.h"

// Manually test insertion, erasure, and retrieval.
TEST(RandomAccessSetTest, InsertAccessErase) {
  dca::util::RandomAccessSet<std::string> set;
  // Set is empty
  EXPECT_FALSE(set.erase("foo"));

  EXPECT_TRUE(set.insert("foo"));
  EXPECT_TRUE(set.insert("bar"));
  EXPECT_FALSE(set.insert("bar"));
  EXPECT_EQ(2, set.size());

  EXPECT_EQ("foo", set.findByIndex(1));  // foo > bar.
  EXPECT_EQ("bar", set.findByIndex(0));

  EXPECT_THROW(set.findByIndex(2), std::logic_error);

  // Set is now empty
  EXPECT_TRUE(set.erase("foo"));
  EXPECT_TRUE(set.erase("bar"));
  EXPECT_EQ(0, set.size());

  // Test insertion after root has been deleted.
  set.insert("baz");
  EXPECT_EQ(1, set.size());
  EXPECT_EQ("baz", set.findByIndex(0));
}

// Perform the test with a number of randomly inserted and removed values.
TEST(RandomAccessSetTest, LinearizeAndRandomAccess) {
  dca::util::RandomAccessSet<int> my_set;
  std::set<int> std_set;

  const int n_insertions = 49;
  const int n_removals = 36;

  // Prepare a shuffled list of unique keys.
  std::vector<int> keys(n_insertions);
  std::iota(keys.begin(), keys.end(), 0);
  std::random_shuffle(keys.begin(), keys.end());

  for (int i = 0; i < n_insertions; ++i) {
    std_set.insert(keys[i]);

    my_set.insert(keys[i]);
    ASSERT_TRUE(my_set.checkConsistency());
  }

  // Remove random keys.
  std::mt19937_64 rng(0);
  for (int i = 0; i < n_removals; ++i) {
    const std::size_t key_idx = std::uniform_int_distribution<std::size_t>(0, keys.size() - 1)(rng);
    const auto key = keys.at(key_idx);
    keys.erase(keys.begin() + key_idx);

    std_set.erase(key);

    my_set.erase(key);
    ASSERT_TRUE(my_set.checkConsistency());
  }

  auto linearized = my_set.linearize();

  ASSERT_EQ(std_set.size(), my_set.size());

  std::size_t idx = 0;
  for (auto it : std_set) {
    EXPECT_EQ(it, linearized[idx]);

    // Test random accessor.
    EXPECT_EQ(it, my_set.findByIndex(idx));

    ++idx;
  }
}

TEST(RandomAccessSetTest, Assignment) {
  dca::util::RandomAccessSet<float> set1{0.5, 3.14, -273.15};

  dca::util::RandomAccessSet<float> set2;
  dca::util::RandomAccessSet<float> set3;

  set2 = set1;
  EXPECT_EQ(set1.linearize(), set2.linearize());

  set3 = std::move(set1);
  EXPECT_EQ(set2.linearize(), set3.linearize());
}
