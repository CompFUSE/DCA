// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Tests call_once_per_loop.hpp

#include "dca/util/call_once_per_loop.hpp"

#include <vector>
#include <gtest/gtest.h>

#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"

void task(uint loop_id, std::vector<int>& data) {
  static dca::util::OncePerLoopFlag flag;
  // to prevent optimizer removing dummy loop
  static volatile int dummy = 0;

  dca::util::callOncePerLoop(flag, loop_id, [&]()
    {
      ++data[loop_id];
      for (int i=0; i<1000; i++) {
        for (int j=0; j<1000; j++) {
          dummy++;
        }
      }
    }
  );
}

TEST(CallOncePerLoopTest, All) {
  const int n_loops(1000);
  std::vector<int> result(n_loops, 0);

  {
    const int n_threads = 8;
    dca::parallel::ThreadPool pool(n_threads);

    for (uint loop_id = 0; loop_id < n_loops; ++loop_id) {
      for (int thread_id = 0; thread_id < n_threads; ++thread_id) {
        pool.enqueue(task, loop_id, std::ref(result));
      }
    }
  }

  for (int elem : result)
    EXPECT_EQ(1, elem);
}
