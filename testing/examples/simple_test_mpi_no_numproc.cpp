#include <iostream>
#include <mpi.h>

#include "gtest/gtest.h"
#include "../common/minimalist_printer.hpp"

TEST(DummyTest1, someTest) {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  EXPECT_EQ(4, size);
}

int main(int argc, char **argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   

  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
