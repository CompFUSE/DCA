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
// This file provides specific tests for the ADIOS2 reader and writer.

#include "dca/io/adios2/adios2_reader.hpp"
#include "dca/io/adios2/adios2_writer.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/testing/minimalist_printer.hpp"

#include <array>
#include <complex>
#include <string>
#include <vector>
#include <typeinfo>  // operator typeid

#include "gtest/gtest.h"

int rank, comm_size;
dca::parallel::MPIConcurrency* concurrency_ptr;

template <typename Scalar>
class ADIOS2ParallelIOTest : public ::testing::Test {};
using TestTypes = ::testing::Types<float, std::complex<double>>;
TYPED_TEST_CASE(ADIOS2ParallelIOTest, TestTypes);

TYPED_TEST(ADIOS2ParallelIOTest, FunctionReadWrite) {
  using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<4>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<2>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  dca::func::function<Scalar, Dmn> f1("parallelFunc");
  size_t dmn_size = 1;
  for (int l = 0; l < f1.signature(); ++l)
    dmn_size *= f1[l];
  int val = rank * dmn_size;
  for (auto& x : f1)
    x = ++val;

  uint64_t start = 0;
  uint64_t end = 0;
  // This returns the linearized bounds of the function for a rank.
  dca::parallel::util::getComputeRange(concurrency_ptr->id(), concurrency_ptr->number_of_processors(),
                                     f1.size(), start, end);
    
  {
    dca::io::ADIOS2Writer writer(concurrency_ptr);
    writer.open_file("test_func_" + typeStr + ".bp", true);

    // Because the caller needs to know if its function is distributed or not we will assume this is so for the API as well.
    // in the future I think something more sophisticated needs to be done and the function will need to know its
    // distribution, but for now we distribute only over the fastest index

    writer.execute(f1, start, end);

    writer.close_file();

    // Read test file.
    /*
    dca::io::ADIOS2Reader reader;
    reader.open_file("test_func.bp");

    dca::func::function<Scalar, Dmn> f2("parallelFunc");

    EXPECT_TRUE(reader.execute(f2));

    for (int i = 0; i < f1.size(); ++i)
      EXPECT_EQ(f1(i), f2(i));

    reader.close_file();
    */
  }
}

int main(int argc, char** argv) {
  int result = 0;

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  rank = concurrency.id();
  comm_size = concurrency.number_of_processors();
  concurrency_ptr = &concurrency;

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
