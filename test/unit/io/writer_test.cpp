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
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "gtest/gtest.h"
#include "dca/util/ignore.hpp"
#include <memory>

template <typename T>
class WriterTest : public ::testing::Test {
};

// dca::io::CSVWriter, dca::io::HDF5Writer, 
using WriterTypes = ::testing::Types<dca::io::HDF5Writer, dca::io::JSONWriter>;
TYPED_TEST_CASE(WriterTest, WriterTypes);

TYPED_TEST(WriterTest, Constructor) {
  TypeParam writer;
  dca::util::ignoreUnused(writer);
}

TYPED_TEST(WriterTest, Unique_Ptr) {
  try {
    TypeParam writer;
    writer.open_file("test_write");
    std::unique_ptr<std::vector<float>> up_vec{new std::vector<float>({1.0,2.0,3.0})};

    writer.execute("test_unique_ptr_to_vector", *up_vec);
    writer.execute("test_unique_ptr_to_vector", std::move(up_vec));
    std::unique_ptr<std::vector<float>> up_vec_null;
    writer.execute("test_unique_ptr_to_vector", std::move(up_vec_null));
  } catch (const std::logic_error& e) {
    ADD_FAILURE() << e.what();
  }
}

