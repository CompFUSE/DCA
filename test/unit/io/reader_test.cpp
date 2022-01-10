// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file tests the readers of the IO library.

#include "dca/io/csv/csv_reader.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include <vector>
#include "gtest/gtest.h"
#include "dca/util/ignore.hpp"

template <typename T>
class ReaderTest : public ::testing::Test {};

using ReaderTypes = ::testing::Types<dca::io::CSVReader, dca::io::HDF5Reader, dca::io::JSONReader
#ifdef DCA_HAVE_ADIOS2
                                     ,dca::io::ADIOS2Reader
#endif
                                     >;
TYPED_TEST_CASE(ReaderTest, ReaderTypes);

TYPED_TEST(ReaderTest, Constructor) {
  TypeParam reader;
  dca::util::ignoreUnused(reader);
}

TEST(JSONReaderTest, Vector) {
  dca::io::JSONReader reader;
  reader.open_file(DCA_SOURCE_DIR "/test/unit/io/vector.json");

  // Simple 3D vector
  std::vector<double> vec_1;
  std::vector<double> vec_1_check{1., 2., 3.};
  reader.execute("simple-vector", vec_1);
  EXPECT_EQ(vec_1_check, vec_1);

  // Vector of two 2D vectors
  std::vector<std::vector<double>> vec_2;
  std::vector<std::vector<double>> vec_2_check{{1.2, 3.4}, {5.6, 7.8}};
  reader.execute("vector-of-vectors", vec_2);
  EXPECT_EQ(vec_2_check, vec_2);

  reader.close_file();
}

TEST(HDF5ReaderTest, Vector) {
  dca::io::HDF5Reader reader;
  reader.open_file(DCA_SOURCE_DIR "/test/unit/io/vector.hdf5");

  // Simple 3D vector
  std::vector<float> vec_1;
  std::vector<float> vec_1_check{1., 2., 3.};
  reader.execute("simple-vector", vec_1);
  EXPECT_EQ(vec_1_check, vec_1);

  // Vector of 3 2D vectors variable length
  std::vector<std::vector<float>> vec_2;
  std::vector<std::vector<float>> vec_2_check{{1.2, 3.4}, {5.6, 7.8, 4.4}, {1.0}};
  reader.execute("vector-of-vectors", vec_2);
  EXPECT_EQ(vec_2_check, vec_2);

  reader.close_file();
}
