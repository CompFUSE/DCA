// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests mpi_collective_sum.hpp.
//
// This test only passes for 8 MPI processes.

#include "dca/parallel/mpi_concurrency/mpi_gather.hpp"

#include <memory>

#include "gtest/gtest.h"

#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/parallel/mpi_concurrency/mpi_gang.hpp"
#include "dca/function/function.hpp"
#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/function/domains/local_domain.hpp"
#include "dca/testing/minimalist_printer.hpp"

std::unique_ptr<dca::parallel::MPIConcurrency> concurrency;

TEST(MPIGatherTest, GatherLocalDmn) {
  using Dmn1 = dca::func::dmn<4>;
  using Dmn2 = dca::func::dmn<10>;
  using LocalDmn = dca::func::LocalDomain<Dmn2>;

  std::vector<int> val1(Dmn1::get_size());
  for (int i = 0; i < val1.size(); ++i)
    val1[i] = i;
  Dmn1::set_elements(val1);

  std::vector<int> val2(Dmn2::get_size());
  for (int i = 0; i < val2.size(); ++i)
    val2[i] = i;
  Dmn2::set_elements(val2);

  dca::parallel::MPIGang gang(*concurrency, 4);

  LocalDmn::initialize(gang);

  dca::func::function<int, dca::func::dmn_variadic<dca::func::dmn_0<Dmn1>, dca::func::dmn_0<LocalDmn>>> local_f;
  dca::func::function<int, dca::func::dmn_variadic<dca::func::dmn_0<Dmn1>, dca::func::dmn_0<Dmn2>>> f;

  for (int i2 = 0; i2 < LocalDmn::get_physical_size(); ++i2)
    for (int i1 = 0; i1 < Dmn1::get_size(); ++i1)
      local_f(i1, i2) = Dmn1::get_elements()[i1] * LocalDmn::get_elements()[i2];


  concurrency->gather(local_f, f, gang);

  for (int i2 = 0; i2 < Dmn2::get_size(); ++i2)
    for (int i1 = 0; i1 < Dmn1::get_size(); ++i1)
      EXPECT_EQ(Dmn1::get_elements()[i1] * Dmn2::get_elements()[i2], f(i1, i2));
}

int main(int argc, char** argv) {
  int result = 0;

  concurrency = std::make_unique<dca::parallel::MPIConcurrency>(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (concurrency->id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
