// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
