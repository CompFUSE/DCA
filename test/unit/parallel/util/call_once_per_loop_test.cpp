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

#include "dca/config/threading.hpp"
#include "dca/parallel/util/call_once_per_loop.hpp"

#include <vector>
#include <gtest/gtest.h>

int task(unsigned int loop_id, std::vector<int>& data) {
  static dca::util::OncePerLoopFlag flag;

  dca::util::callOncePerLoop(flag, loop_id, [&]() {
    ++data[loop_id];
#ifdef DCA_HAVE_HPX
    hpx::this_thread::sleep_for(std::chrono::microseconds(100));
#else
    std::this_thread::sleep_for(std::chrono::microseconds(100));
#endif
  });
  return 0;
}

TEST(CallOncePerLoopTest, All) {
  const int n_loops(1000);
  std::vector<int> result(n_loops, 0);

  {
    const int n_threads = 8;
#ifdef DCA_HAVE_HPX
#else
    dca::parallel::ThreadPool pool(n_threads);
#endif
    for (unsigned int loop_id = 0; loop_id < n_loops; ++loop_id) {
    std::vector<dca::parallel::thread_traits::future_type<int>> futures;
      for (int thread_id = 0; thread_id < n_threads; ++thread_id) {
#ifdef DCA_HAVE_HPX
        futures.push_back(hpx::async(task, loop_id, std::ref(result)));
#else
	pool.enqueue(task, loop_id, std::ref(result));
#endif
      }
#ifdef DCA_HAVE_HPX
        hpx::wait_all(futures);
#endif
    }

  }

  for (int elem : result)
    EXPECT_EQ(1, elem);
}
