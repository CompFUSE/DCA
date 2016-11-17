// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests print_time.hpp.

#include "dca/util/print_time.hpp"
#include <iostream>
#include "gtest/gtest.h"

TEST(PrintTimeTest, std_chrono_system_clock) {
  std::cout << dca::util::print_time() << std::endl;
}
