// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests timer.hpp.

#include "dca/util/timer.hpp"
#include <thread>
#include "gtest/gtest.h"

TEST(TimerTest, PrintingThread) {
  dca::util::Timer timer("sleep-for-1ms");

  using namespace std::chrono_literals;
  std::this_thread::sleep_for(1ms);
}

TEST(TimerTest, SilentThread) {
  dca::util::Timer timer("sleep-for-1ms", false);
}
