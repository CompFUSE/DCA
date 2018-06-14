// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// For the moment this file just tests the construction of the different types of writers.

#include "gtest/gtest.h"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/io/csv/csv_writer.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/util/ignore.hpp"

#include <memory>

template <typename T>
class WriterTest : public ::testing::Test {};

using WriterTypes =
    ::testing::Types<dca::io::HDF5Writer, dca::io::JSONWriter>;  //, dca::io::CSVWriter>;
TYPED_TEST_CASE(WriterTest, WriterTypes);

TYPED_TEST(WriterTest, Constructor) {
  TypeParam writer;
  dca::util::ignoreUnused(writer);
}

TYPED_TEST(WriterTest, Unique_Ptr) {
  TypeParam writer;
  std::string test_name(::testing::UnitTest::GetInstance()->current_test_info()->name());
  std::string type_string(::testing::UnitTest::GetInstance()->current_test_info()->type_param());

  writer.open_file("/tmp/test_file_" + type_string + "_" + test_name);
  writer.open_group("writer_test");

  std::unique_ptr<std::vector<float>> up_vec{new std::vector<float>({1.0, 2.0, 3.0, 4.0})};

  writer.execute("test_unique_ptr_to_vector", up_vec);
  std::unique_ptr<std::vector<float>> up_vec_null;
  writer.execute("test_null_unique_ptr_to_vector", up_vec_null);

  using test_domain_0a = dca::func::dmn_0<dca::func::dmn<1, double>>;
  using TestFunc = dca::func::function<double, test_domain_0a>;

  std::unique_ptr<TestFunc> up_test_func{new TestFunc("Test_Func")};
  up_test_func->reset();
  writer.execute(up_test_func);
  std::unique_ptr<TestFunc> test_function_null = nullptr;
  writer.execute(test_function_null);

  writer.close_group();
  writer.close_file();
}
