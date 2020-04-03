// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides a map with O(log n) insertion, removal and random access (access of the value relative
// to the i-th lowest key. Useful for a random selection in an ordered list with variable size.

#include "dca/util/random_access_map.hpp"

#include <map>
#include <random>
#include <string>

#include "gtest/gtest.h"

// Manually test insertion, erasure, and retrieval.
TEST(RandomAccessMapTest, InsertFindErase) {
  dca::util::RandomAccessMap<std::string, int> map;
  // Map is empty
  EXPECT_THROW(map.erase("foo"), std::logic_error);

  map.insert("foo", 2);
  map.insert("bar", 1);
  EXPECT_EQ(2, map.size());

  EXPECT_EQ(2, map.find("foo"));
  EXPECT_EQ(1, map.find("bar"));
  EXPECT_EQ(2, map[1]);  // foo > bar.
  EXPECT_EQ(1, map[0]);
  EXPECT_THROW(map.find("baz"), std::logic_error);
  EXPECT_THROW(map[2], std::logic_error);

  // Map is now empty
  map.erase("foo");
  map.erase("bar");
  EXPECT_EQ(0, map.size());

  // Test insertion after root has been deleted.
  map.insert("baz", 3);
  EXPECT_EQ(1, map.size());
  EXPECT_EQ(3, map.find("baz"));
  EXPECT_EQ(3, map[0]);
}

// Perform the test with a number of randomly inserted and removed values.
TEST(RandomAccessMapTest, LinearizeAndRandomAccess) {
  dca::util::RandomAccessMap<int, int> my_map;
  std::map<int, int> std_map;

  const int n_insertions = 100;
  const int n_removals = 75;

  // Prepare a shuffled list of unique keys.
  std::vector<int> keys(n_insertions);
  std::iota(keys.begin(), keys.end(), 0);
  std::random_shuffle(keys.begin(), keys.end());

  for (int i = 0; i < n_insertions; ++i) {
    const int val = i;
    std_map[keys[i]] = val;

    my_map.insert(keys[i], val);
    ASSERT_TRUE(my_map.checkConsistency());
  }

  // Remove random keys.
  std::mt19937_64 rng(0);
  for (int i = 0; i < n_removals; ++i) {
    const std::size_t key_idx = std::uniform_int_distribution<std::size_t>(0, keys.size() - 1)(rng);
    const auto key = keys.at(key_idx);
    keys.erase(keys.begin() + key_idx);

    std_map.erase(key);

    my_map.erase(key);
    ASSERT_TRUE(my_map.checkConsistency());
  }

  auto linearized = my_map.linearize();

  ASSERT_EQ(std_map.size(), my_map.size());

  std::size_t idx = 0;
  for (auto it : std_map) {
    EXPECT_EQ(it.first, linearized[idx].first);
    EXPECT_EQ(it.second, linearized[idx].second);

    // Test random accessor.
    EXPECT_EQ(it.second, my_map[idx]);

    ++idx;
  }
}

TEST(RandomAccessMapTest, Assignment) {
  dca::util::RandomAccessMap<int, double> map1{{1, 0.5}, {-1, 3.14}, {42, -273.15}};
  dca::util::RandomAccessMap<int, double> map2;
  dca::util::RandomAccessMap<int, double> map3;

  map2 = map1;
  EXPECT_EQ(map1.linearize(), map2.linearize());

  map3 = std::move(map1);
  EXPECT_EQ(map2.linearize(), map3.linearize());
}
