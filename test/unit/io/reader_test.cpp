// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// For the moment this file just tests the construction of the different types of readers.

#include "dca/io/csv/csv_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "gtest/gtest.h"
#include "dca/util/ignore.hpp"

template <typename T>
class ReaderTest : public ::testing::Test {};

using ReaderTypes = ::testing::Types<dca::io::CSVReader, dca::io::JSONReader>;
TYPED_TEST_CASE(ReaderTest, ReaderTypes);

TYPED_TEST(ReaderTest, Constructor) {
  TypeParam reader;
  dca::util::ignoreUnused(reader);
}
