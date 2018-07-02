// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests stdthread.hpp.

#include "dca/parallel/stdthread/stdthread.hpp"

#include "gtest/gtest.h"

using Threading = dca::parallel::stdthread;

void routine(const int id, const int num_threads, std::vector<int>& vec) {
  EXPECT_EQ(4, num_threads);

  vec[id] += id;
}

TEST(StdthreadTest, Execute) {
  Threading threading;

  const int num_threads = 4;
  std::vector<int> vec{0, 10, 20, 30};
  std::vector<int> vec_check{0, 11, 22, 33};

  threading.execute(num_threads, routine, std::ref(vec));

  EXPECT_EQ(vec_check, vec);
}
