// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the HDF5 reader.

#include "dca/io/hdf5/hdf5_reader.hpp"
#include <string>
#include "gtest/gtest.h"
#include "dca/io/hdf5/hdf5_writer.hpp"

TEST(HDF5ReaderWriterTest, ReaderDestructorCleanUp) {
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

TEST(HDF5ReaderWriterTest, WriterDestructorCleanUp) {
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

TEST(HDF5ReaderWriterTest, VectorReadWrite) {
  const std::string object_name = "a_vector";
  const std::string file_name = "hdf5_reader_vector_test.hdf5";
  const std::vector<std::complex<double>> a_vector{
      std::complex<double>(1., 0.), std::complex<double>(0., 1.), std::complex<double>(23.4, -1.5)};

  // Create test file.
  dca::io::HDF5Writer writer;
  writer.open_file(file_name);
  writer.execute(object_name, a_vector);
  writer.close_file();

  // Read test file.
  dca::io::HDF5Reader reader;
  std::vector<std::complex<double>> vector_read;
  reader.open_file(file_name);
  reader.execute(object_name, vector_read);

  ASSERT_EQ(a_vector.size(), vector_read.size());
  for (int i = 0; i < a_vector.size(); ++i) {
    EXPECT_DOUBLE_EQ(std::real(a_vector[i]), std::real(vector_read[i]));
    EXPECT_DOUBLE_EQ(std::imag(a_vector[i]), std::imag(vector_read[i]));
  }

  reader.close_file();
}
