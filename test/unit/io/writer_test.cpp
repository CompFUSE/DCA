// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// For the moment this file just tests the construction of the different types of writers.

#include "dca/io/csv/csv_writer.hpp"
#include "gtest/gtest.h"
#include "dca/util/ignore.hpp"

template <typename T>
class WriterTest : public ::testing::Test {};

using WriterTypes = ::testing::Types<dca::io::CSVWriter>;
TYPED_TEST_CASE(WriterTest, WriterTypes);

TYPED_TEST(WriterTest, Constructor) {
  TypeParam writer;
  dca::util::ignoreUnused(writer);
}
