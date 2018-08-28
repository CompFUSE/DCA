// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests mpi_packing.hpp.
// It is run with only 1 MPI process since packing does not involve any communication.

#include "dca/parallel/mpi_concurrency/mpi_packing.hpp"

#include <cstring>  // std::memcpy

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/testing/minimalist_printer.hpp"

TEST(MPIPackingTest, PackAndUnpackFunction) {
  dca::parallel::MPIPacking packing;

  using TestDomain = dca::func::dmn_0<dca::func::dmn<4, int>>;

  dca::func::function<double, TestDomain> f("test-function");
  for (std::size_t i = 0; i < f.size(); ++i)
    f(i) = i + 3.14;

  //
  // Non-const function
  //
  dca::func::function<double, TestDomain> f_non_const(f, "non-const-function");

  const auto buffer_size_non_const = packing.get_buffer_size(f_non_const);
  char* buffer_non_const = new char[buffer_size_non_const];
  int offset_non_const = 0;

  packing.pack(buffer_non_const, buffer_size_non_const, offset_non_const, f_non_const);

  EXPECT_EQ(buffer_size_non_const, offset_non_const);
  for (std::size_t i = 0; i < f.size(); ++i)
    EXPECT_EQ(f(i), f_non_const(i));

  dca::func::function<double, TestDomain> f_unpacked("unpacked");

  // Use a copy of the buffer in unpack to detect if pack has written out of bounds.
  char* buffer_non_const_copy = new char[buffer_size_non_const];
  std::memcpy(buffer_non_const_copy, buffer_non_const, buffer_size_non_const);

  int offset_unpack = 0;
  packing.unpack(buffer_non_const_copy, buffer_size_non_const, offset_unpack, f_unpacked);

  EXPECT_EQ(buffer_size_non_const, offset_unpack);
  for (std::size_t i = 0; i < f.size(); ++i)
    EXPECT_EQ(f(i), f_unpacked(i));

  delete[] buffer_non_const;
  delete[] buffer_non_const_copy;

  //
  // Const function
  //
  const dca::func::function<double, TestDomain> f_const(f, "const-function");

  const auto buffer_size_const = packing.get_buffer_size(f_const);
  char* buffer_const = new char[buffer_size_const];
  int offset_const = 0;

  packing.pack(buffer_const, buffer_size_const, offset_const, f_const);

  EXPECT_EQ(buffer_size_const, offset_const);

  char* buffer_const_copy = new char[buffer_size_const];
  std::memcpy(buffer_const_copy, buffer_const, buffer_size_const);

  offset_unpack = 0;
  packing.unpack(buffer_const_copy, buffer_size_const, offset_unpack, f_unpacked);

  EXPECT_EQ(buffer_size_const, offset_unpack);
  for (std::size_t i = 0; i < f_const.size(); ++i)
    EXPECT_EQ(f_const(i), f_unpacked(i));

  delete[] buffer_const;
  delete[] buffer_const_copy;
}

int main(int argc, char** argv) {
  int result = 0;

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
