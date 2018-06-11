// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author:  Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests mpi_collective_gather.hpp.

#include "dca/parallel/mpi_concurrency/mpi_collective_gather.hpp"

#include <complex>
#include <memory>

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/testing/minimalist_printer.hpp"

class MPICollectiveGatherTest : public ::testing::Test {
protected:
  MPICollectiveGatherTest() {
    grouping_.reset(new dca::parallel::MPIProcessorGrouping);
    size_ = grouping_->get_Nr_threads();
    rank_ = grouping_->get_id();
  }

  int size_;
  int rank_;

  std::unique_ptr<const dca::parallel::MPIProcessorGrouping> grouping_;
};

TEST_F(MPICollectiveGatherTest, GatherVector) {
  auto some_work = [](const int id) { return std::complex<float>(id, 0); };

  std::vector<std::complex<float>> v(10);
  const auto bounds = dca::parallel::util::getBounds(rank_, size_, std::make_pair(0, v.size()));

  for (int i = bounds.first; i < bounds.second; ++i)
    v[i] = some_work(i);

  dca::parallel::MPICollectiveGather gather_interface(grouping_);
  gather_interface.gather(v);

  for (int i = 0; i < v.size(); ++i)
    EXPECT_EQ(some_work(i), v[i]);
}

TEST_F(MPICollectiveGatherTest, GatherFunction) {
  auto some_work = [](const int id) { return id * id; };

  dca::func::function<int, dca::func::dmn_variadic<dca::func::dmn_0<dca::func::dmn<4>>,
                                                   dca::func::dmn_0<dca::func::dmn<7>>>>
      f;

  const auto bounds = dca::parallel::util::getBounds(rank_, size_, std::make_pair(0, f.size()));

  for (int i = bounds.first; i < bounds.second; ++i)
    f(i) = some_work(i);

  dca::parallel::MPICollectiveGather gather_interface(grouping_);
  gather_interface.gather(f);

  for (int i = 0; i < f.size(); ++i)
    EXPECT_EQ(some_work(i), f(i));
}

int main(int argc, char** argv) {
  int result = 0;
  int rank;

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
