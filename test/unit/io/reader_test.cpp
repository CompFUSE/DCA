// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
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
#include <functional>
#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_reader.hpp"
#endif
#include <vector>
#include "gtest/gtest.h"
#include "dca/util/ignore.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"

template <typename T>
class ReaderTest : public ::testing::Test {};


#ifdef DCA_HAVE_MPI
using Concurrency = dca::parallel::MPIConcurrency;
#else
using Concurrency = dca::parallel::NoConcurrency;
#endif




using ReaderTypes = ::testing::Types<dca::io::CSVReader, dca::io::HDF5Reader, dca::io::JSONReader
#ifdef DCA_HAVE_ADIOS2
                                     ,dca::io::ADIOS2Reader<Concurrency>
#endif
                                     >;

TYPED_TEST_CASE(ReaderTest, ReaderTypes);

std::unique_ptr<Concurrency> concurrency_ptr;

TEST(JSONReaderTest, Vector) {
  dca::io::JSONReader reader;
  reader.open_file(DCA_SOURCE_DIR "/test/unit/io/vector.json");

  // Simple 3D vector
  std::vector<double> vec_1;
  std::vector<double> vec_1_check{1., 2., 3.};
  reader.execute("simple-vector", vec_1);
  EXPECT_EQ(vec_1_check, vec_1);

  // Vector of two vectors of variable length
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
  bool result = reader.execute("simple-vector", vec_1);
  EXPECT_TRUE(result);
  EXPECT_EQ(vec_1_check, vec_1);

  // Vector of 3 vectors variable length
  std::vector<std::vector<float>> vec_2;
  std::vector<std::vector<float>> vec_2_check{{1.2, 3.4}, {5.6, 7.8, 4.4}, {1.0}};
  result = reader.execute("vector-of-vectors", vec_2);
  EXPECT_TRUE(result);
  EXPECT_EQ(vec_2_check, vec_2);

  reader.close_file();
}

#ifdef DCA_HAVE_ADIOS2
// more extensive testing in adios2_reader_writer_test
TEST(ADIOS2ReaderTest, Vector) {
  dca::io::ADIOS2Reader reader(*concurrency_ptr, true);
  reader.open_file(DCA_SOURCE_DIR "/test/unit/io/vector.bp");

  // Simple 3D vector
  std::vector<double> vec_1;
  std::vector<double> vec_1_check{1., 2., 3.};
  bool found = reader.execute("simple-vector", vec_1);
  EXPECT_TRUE(found);
  EXPECT_EQ(vec_1_check, vec_1);

  // vector of 2 vector of constant length (ADIOS can't handle variable length collections)
  // it can handle variables whose shape changes step to step.
  std::vector<std::vector<float>> vec_2;
  std::vector<std::vector<float>> vec_2_check{{1.2, 3.4}, {5.6, 7.8, 4.4}, {1.0}};
  reader.execute("vector-of-vectors", vec_2);
  EXPECT_EQ(vec_2_check, vec_2);

  reader.close_file();
}
#endif

int main(int argc, char** argv) {
  concurrency_ptr = std::make_unique<Concurrency>(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
