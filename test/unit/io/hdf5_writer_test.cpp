// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides specific tests for the HDF5 writer.

#include "dca/io/hdf5/hdf5_writer.hpp"
#include <string>
#include "gtest/gtest.h"

TEST(HDF5WriterTest, DestructorCleanUp) {
  std::string test_file_name = "hdf5_writer_test.hdf5";
  std::string group_name_1 = "integers";
  std::string group_name_2 = "magic-numbers";
  std::string object_name = "forty-two";

  dca::io::HDF5Writer writer;
  const int i = 42;

  writer.open_file(test_file_name);
  writer.open_group(group_name_1);
  writer.open_group(group_name_2);
  writer.execute(object_name, i);

  // Destructor closes groups and file.
  // writer.close_group();
  // writer.close_group();
  // writer.close_file();
}
