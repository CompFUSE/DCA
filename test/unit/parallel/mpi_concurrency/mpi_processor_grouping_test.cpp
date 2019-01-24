// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests mpi_processor_grouping.hpp.
// It is run with 4 MPI processes.

#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "gtest/gtest.h"
#include "dca/testing/minimalist_printer.hpp"

// Global variable that can be accessed in the test.
int rank;

TEST(MPIProcessorGroupingTest, All) {
  dca::parallel::MPIProcessorGrouping grouping;

  EXPECT_EQ(rank, grouping.get_id());
  EXPECT_EQ(4, grouping.get_size());
  EXPECT_EQ(0, grouping.first());
  EXPECT_EQ(3, grouping.last());
}

TEST(MPIProcessorGroupingTest, FaultyProcess) {
  // Simulate a faulty driver at ranks 2 and 0.
  auto mock_check = [] { return rank != 2 && rank != 0; };

  dca::parallel::MPIProcessorGrouping grouping(mock_check);

  if (rank == 2 || rank == 0) {
    EXPECT_FALSE(grouping.isValid());
  }
  else {
    EXPECT_TRUE(grouping.isValid());

    if (rank == 1)
      EXPECT_EQ(0, grouping.get_id());
    else if (rank == 3)
      EXPECT_EQ(1, grouping.get_id());

    EXPECT_EQ(2, grouping.get_size());
    EXPECT_EQ(0, grouping.first());
    EXPECT_EQ(1, grouping.last());
  }
}

int main(int argc, char** argv) {
  int result = 0;

  MPI_Init(&argc, &argv);
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
