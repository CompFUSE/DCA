// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests serial_processor_grouping.hpp.

#include "dca/parallel/no_concurrency/serial_processor_grouping.hpp"
#include "gtest/gtest.h"

TEST(SerialProcessorGroupingTest, All) {
  dca::parallel::SerialProcessorGrouping grouping;

  EXPECT_EQ(0, grouping.get_id());
  EXPECT_EQ(1, grouping.get_Nr_threads());
  EXPECT_EQ(0, grouping.first());
  EXPECT_EQ(0, grouping.last());
}
