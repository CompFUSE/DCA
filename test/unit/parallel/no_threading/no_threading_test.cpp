// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests no_threading.hpp.

#include "dca/parallel/no_threading/no_threading.hpp"
#include "gtest/gtest.h"

void* start_routine(void* arg) {
  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(arg);

  const int id = static_cast<int>(data_ptr->id);
  const int num_threads = static_cast<int>(data_ptr->num_threads);
  std::vector<int>* vec_ptr = static_cast<std::vector<int>*>(data_ptr->arg);

  EXPECT_EQ(4, num_threads);

  vec_ptr->operator[](id) += id;

  return 0;
}

TEST(NoThreadingTest, Execute) {
  dca::parallel::NoThreading threading;

  const int num_threads = 4;
  std::vector<int> vec{0, 10, 20, 30};
  std::vector<int> vec_check{0, 11, 22, 33};

  threading.execute(num_threads, start_routine, static_cast<void*>(&vec));

  EXPECT_EQ(vec_check, vec);
}

TEST(NoThreadingTest, OstreamOperator) {
  dca::parallel::NoThreading threading;
  const int num_threads = 4;
  std::vector<int> vec{0, 10, 20, 30};
  threading.execute(num_threads, start_routine, static_cast<void*>(&vec));
  std::string eout("\nthreading type:NoThreading\nnumber of threads:4");
  EXPECT_EQ(eout, ::testing::PrintToString(threading));
}
