// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file the collective sum of a large function object.
// Warning: 8 GB of memory will be allocated.

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

// Note: the minimalist printer is not used for compatibility with the Summit cluster at the moment of writing.

using dca::func::dmn;
using dca::func::dmn_0;
using dca::func::dmn_variadic;
using dca::func::function;

TEST(MPICollectiveSumTest, SumScalar) {
  dca::parallel::MPIConcurrency concurrency(0, nullptr);
  const int rank = concurrency.id();
  const int comm_size = concurrency.number_of_processors();

  using Dmn1 = dmn_0<dmn<(1 << 16), int>>;
  using Dmn2 = dmn_0<dmn<2 + (1 << 15), char>>;

  dca::func::function<std::uint8_t, dmn_variadic<Dmn1, Dmn2>> large_f;
  ASSERT_EQ(static_cast<unsigned long int>(Dmn1::dmn_size()) * Dmn2::dmn_size(), large_f.size());
  ASSERT_GT(large_f.size(), std::numeric_limits<int>::max());

  const int rank_sum = (comm_size * (comm_size + 1)) / 2;  // considers the +1 in (rank + 1).
  ASSERT_LE(2 * rank_sum, std::numeric_limits<std::uint8_t>::max());

  for (std::size_t i = 0; i < large_f.size(); i+=2) {
    large_f(i) = (rank + 1) * 1;
    large_f(i + 1) = (rank + 1) * 2;
  }

  concurrency.sum(large_f);

  for (std::size_t i = 0; i < large_f.size(); i+=2) {
    EXPECT_EQ(1 * rank_sum, large_f(i));
    EXPECT_EQ(2 * rank_sum, large_f(i + 1));
  }
}
