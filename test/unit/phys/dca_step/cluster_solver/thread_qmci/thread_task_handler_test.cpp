// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the thread task handler.

#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include <numeric>
#include <string>
#include <vector>
#include "gtest/gtest.h"

namespace dca {
namespace testing {
// dca::testing::

void constructorTestBody(const int num_walkers, const int num_accumulators,
                         const std::vector<std::string>& expected_thread_tasks) {
  dca::phys::solver::ThreadTaskHandler handler(num_walkers, num_accumulators);

  EXPECT_EQ(expected_thread_tasks, handler.getThreadTasksVec());

  for (int id = 0; id < expected_thread_tasks.size(); ++id) {
    EXPECT_EQ(expected_thread_tasks[id], handler.getTask(id));
  }
}

void walkerIDToRngIndexTestBody(const int num_walkers, const int num_accumulators) {
  dca::phys::solver::ThreadTaskHandler handler(num_walkers, num_accumulators);

  std::vector<int> rng_indices;
  std::vector<int> expected(num_walkers);
  std::iota(expected.begin(), expected.end(), 0);

  for (int id = 0; id < num_walkers + num_accumulators; ++id) {
    if (handler.getTask(id) == "walker")
      rng_indices.push_back(handler.walkerIDToRngIndex(id));
  }
  EXPECT_EQ(expected, rng_indices);
}

}  // testing
}  // dca

TEST(ThreadTaskHandlerTest, Constructor) {
  // 0 walkers, 0 accumulators
  std::vector<std::string> expected = {};
  dca::testing::constructorTestBody(0, 0, expected);

  // 1 walker, 0 accumulators
  expected = {"walker"};
  dca::testing::constructorTestBody(1, 0, expected);

  // 0 walkers, 1 accumulator
  expected = {"accumulator"};
  dca::testing::constructorTestBody(0, 1, expected);

  // 1 walker, 1 accumulator
  expected = {"walker", "accumulator"};
  dca::testing::constructorTestBody(1, 1, expected);

  // 4 walkers, 2 accumulators
  expected = {"walker", "accumulator", "walker", "accumulator", "walker", "walker"};
  dca::testing::constructorTestBody(4, 2, expected);

  // 2 walkers, 4 accumulators
  expected = {"walker", "accumulator", "walker", "accumulator", "accumulator", "accumulator"};
  dca::testing::constructorTestBody(2, 4, expected);
}

TEST(ThreadTaskHandlertest, walkerIDToRngIndex) {
  // 0 walkers, 0 accumulators
  dca::testing::walkerIDToRngIndexTestBody(0, 0);

  // 1 walker, 0 accumulators
  dca::testing::walkerIDToRngIndexTestBody(1, 0);

  // 0 walkers, 1 accumulator
  dca::testing::walkerIDToRngIndexTestBody(0, 1);

  // 1 walker, 1 accumulator
  dca::testing::walkerIDToRngIndexTestBody(1, 1);

  // 4 walkers, 2 accumulators
  dca::testing::walkerIDToRngIndexTestBody(4, 2);

  // 2 walkers, 4 accumulators
  dca::testing::walkerIDToRngIndexTestBody(2, 4);
}

#ifndef NDEBUG
TEST(ThreadTaskHandlerDeathTest, walkerIDToRngIndex) {
  dca::phys::solver::ThreadTaskHandler handler(4, 2);  // = w, a, w, a, w, w

  // Call with thread id that is out of bound.
  EXPECT_DEATH(handler.walkerIDToRngIndex(-1),
               "walker_id >= 0 && walker_id < thread_tasks_.size()");
  EXPECT_DEATH(handler.walkerIDToRngIndex(6), "walker_id >= 0 && walker_id < thread_tasks_.size()");

  // Call with thread id that belongs to an accumulator.
  EXPECT_DEATH(handler.walkerIDToRngIndex(3), "thread_tasks_.walker_id. == .walker.");
}
#endif  // NDEBUG
