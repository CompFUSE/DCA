// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests timer.hpp.

#include "dca/util/timer.hpp"
#include <thread>
#include "gtest/gtest.h"

TEST(TimerTest, SimpleTest) {
  dca::util::Timer timer("sleep-for-1ms");

  using namespace std::chrono_literals;
  std::this_thread::sleep_for(1ms);
}
