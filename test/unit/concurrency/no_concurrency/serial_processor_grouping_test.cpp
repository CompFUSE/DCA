// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests serial_processor_grouping.hpp.

#include "dca/concurrency/no_concurrency/serial_processor_grouping.hpp"
#include "gtest/gtest.h"

TEST(SerialProcessorGroupingTest, All) {
  dca::concurrency::SerialProcessorGrouping grouping;

  EXPECT_EQ(0, grouping.get_id());
  EXPECT_EQ(1, grouping.get_Nr_threads());
  EXPECT_EQ(0, grouping.first());
  EXPECT_EQ(0, grouping.last());
}
