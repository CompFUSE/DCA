// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides specific tests for the HDF5 reader.

#include "dca/io/hdf5/hdf5_reader.hpp"
#include <string>
#include "gtest/gtest.h"
#include "dca/io/hdf5/hdf5_writer.hpp"

TEST(HDF5ReaderTest, DestructorCleanUp) {
  std::string test_file_name = "hdf5_reader_test.hdf5";
  std::string group_name = "magic-numbers";
  std::string object_name = "forty-two";

  // Create test file.
  dca::io::HDF5Writer writer;
  const int i = 42;

  writer.open_file(test_file_name);
  writer.open_group(group_name);
  writer.execute(object_name, i);
  writer.close_group();
  writer.close_file();

  // Read test file.
  dca::io::HDF5Reader reader;
  int j;

  reader.open_file(test_file_name);
  reader.open_group(group_name);
  reader.execute(object_name, j);

  EXPECT_EQ(i, j);

  // Destructor closes file.
  // reader.close_file();
}
