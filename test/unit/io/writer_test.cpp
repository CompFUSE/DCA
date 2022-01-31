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
// For the moment this file just tests the construction of the different types of writers.

#include "gtest/gtest.h"
#include <memory>
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/io/csv/csv_writer.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/util/ignore.hpp"
#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
dca::parallel::MPIConcurrency* mpi_concurrency_ptr = nullptr;
#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_writer.hpp"
#include "dca/io/adios2/adios2_global.hpp"
adios2::ADIOS* adios_ptr = nullptr;
#endif
#endif
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
dca::parallel::NoConcurrency* concurrency_ptr = nullptr;

template <typename T>
class WriterTest : public ::testing::Test {};

using Concurrency = dca::parallel::NoConcurrency;
using MPIConcurrency = dca::parallel::MPIConcurrency;

using WriterTypes = ::testing::Types<dca::io::HDF5Writer, dca::io::JSONWriter
#ifdef DCA_HAVE_ADIOS2
                                     ,
                                     dca::io::ADIOS2Writer<MPIConcurrency>
#endif
                                     >;  //, dca::io::CSVWriter>;
TYPED_TEST_CASE(WriterTest, WriterTypes);

TYPED_TEST(WriterTest, Constructor) {
  std::unique_ptr<TypeParam> writer_ptr;
#ifdef DCA_HAVE_ADIOS2
  if constexpr (std::is_same<TypeParam, dca::io::ADIOS2Writer<MPIConcurrency>>::value)
    writer_ptr = std::make_unique<TypeParam>(*adios_ptr, mpi_concurrency_ptr);
  else
#endif
    writer_ptr = std::make_unique<TypeParam>();
  dca::util::ignoreUnused(writer_ptr);
}

TYPED_TEST(WriterTest, Unique_Ptr) {
  std::unique_ptr<TypeParam> writer_ptr;
#ifdef DCA_HAVE_ADIOS2
  if constexpr (std::is_same<TypeParam, dca::io::ADIOS2Writer<MPIConcurrency>>::value)
    writer_ptr = std::make_unique<TypeParam>(*adios_ptr, mpi_concurrency_ptr);
  else
#endif
    writer_ptr = std::make_unique<TypeParam>();

  TypeParam& writer = *writer_ptr;
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

int main(int argc, char** argv) {
  dca::parallel::NoConcurrency concurrency(argc, argv);

  concurrency_ptr = &concurrency;

#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency mpi_concurrency(argc, argv);
  mpi_concurrency_ptr = &mpi_concurrency;
#ifdef DCA_HAVE_ADIOS2
  adios2::ADIOS adios("", mpi_concurrency_ptr->get());
  adios_ptr = &adios;
#endif
#endif

  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
