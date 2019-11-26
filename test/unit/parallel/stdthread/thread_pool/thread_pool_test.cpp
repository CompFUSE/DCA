// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests thread_pool.hpp.

#include "dca/config/threading.hpp"

#include <numeric>

#include "gtest/gtest.h"

TEST(ThreadPoolTest, Enqueue) {
  const int n_items = 9;
  const int n_threads = 4;
  std::vector<int> input(n_items);
  std::vector<int> output(n_items, 0);
  std::iota(input.begin(), input.end(), 0);

  auto workload = [](const std::size_t id, const std::vector<int>& inp, std::vector<int>& out) {
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    out[id] = inp[id] * inp[id];
  };

  const int n_immediate_checks = 4;
  {
    dca::parallel::ThreadPool pool(n_threads);
    EXPECT_EQ(n_threads, pool.size());

    auto task = std::bind(workload, std::placeholders::_1, std::ref(input), std::ref(output));
    std::vector<dca::parallel::thread_traits::future_type<void>> futures;
    for (std::size_t id = 0; id < n_items; ++id)
      futures.emplace_back(pool.enqueue(task, id));

    // Check the synchronization with futures.
    for (std::size_t id = 0; id < n_immediate_checks; ++id) {
      futures[id].wait();
      EXPECT_EQ(input[id] * input[id], output[id]);
    }
  }

  // Check that the other tasks finished before the pool is destroyed.
  for (std::size_t id = n_immediate_checks; id < n_items; ++id)
    EXPECT_EQ(input[id] * input[id], output[id]);
}

TEST(ThreadPoolTest, Enlarge) {
  dca::parallel::ThreadPool& pool = dca::parallel::ThreadPool::get_instance();
  EXPECT_EQ(std::size_t(0), pool.size());

  pool.enlarge(3);
  EXPECT_EQ(std::size_t(3), pool.size());

  pool.enlarge(1);
  EXPECT_EQ(std::size_t(3), pool.size());

  // Dispatch some work to test if queue enlarging breaks running tasks.
  auto workload = [](int id) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    return id;
  };
  std::vector<dca::parallel::thread_traits::future_type<int>> futures;

  for (std::size_t i = 0; i < 5; ++i)
    futures.emplace_back(pool.enqueue(workload, i));

  pool.enlarge(5);
  EXPECT_EQ(std::size_t(5), pool.size());

  for (int i = 5; i < 12; ++i)
    futures.emplace_back(pool.enqueue(workload, i));

  for (std::size_t i = 0; i < 12; ++i)
    EXPECT_EQ(i, futures[i].get());
}
