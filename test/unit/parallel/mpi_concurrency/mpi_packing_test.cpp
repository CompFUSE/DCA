// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests mpi_packing.hpp.
// It is run with only 1 MPI process since packing does not involve any communication.

#include "dca/parallel/mpi_concurrency/mpi_packing.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/testing/minimalist_printer.hpp"

TEST(MPIPackingTest, GetBufferSizeFunction) {
  dca::parallel::MPIProcessorGrouping grouping;
  grouping.set();

  dca::parallel::MPIPacking packing(grouping);

  using TestDomain = dca::func::dmn_0<dca::func::dmn<4, int>>;

  // Non-const function
  dca::func::function<double, TestDomain> f_non_const("non-const-function");
  // 1 int (number of elements) + 4 doubles (elements) = 36 bytes.
  EXPECT_EQ(36, packing.get_buffer_size(f_non_const));

  // Const function
  const dca::func::function<double, TestDomain> f_const("const-function");
  EXPECT_EQ(36, packing.get_buffer_size(f_const));
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
