// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific
// publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

#ifndef DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCKRNG_HPP
#define DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCKRNG_HPP

#include <iostream>
#include <pthread.h>
#include <vector>
#include "gtest/gtest.h"

class MockRng {
public:
  MockRng() {}
  ~MockRng() {}
  void init_from_id(int id, int max_id);

  static void initialize() { pthread_mutex_init(&rng_lock, NULL); }
  static void finalize() { pthread_mutex_destroy(&rng_lock); }
  static void reset() {
    ids.clear();
    max_ids.clear();
  }
  static void checkSeeds(int n_walk);

private:
  static pthread_mutex_t rng_lock;
  static std::vector<int> ids, max_ids;
public:
  int id;
};
pthread_mutex_t MockRng::rng_lock;
std::vector<int> MockRng::ids, MockRng::max_ids;

void MockRng::init_from_id(int walker_id, int max_id) {
  id = walker_id;
  pthread_mutex_lock(&rng_lock);
  ids.push_back(id);
  max_ids.push_back(max_id);
  pthread_mutex_unlock(&rng_lock);
}

void MockRng::checkSeeds(int n_walk) {
  assert(ids.size() == max_ids.size());
  ASSERT_EQ(ids.size(), n_walk);
  std::cout << "Id called, with max " << max_ids[0] << ":\n";
  for (int i = 0; i < ids.size(); i++)
    std::cout << ids[i] << "  ";
  std::cout << std::endl;
  for (int i = 0; i < ids.size(); i++)
    for (int j = i + 1; j < ids.size(); j++)
      ASSERT_NE(ids[i], ids[j]);
  for (int i = 1; i < max_ids.size(); i++)
    ASSERT_EQ(max_ids[i], max_ids[i - 1]);
}

#endif
