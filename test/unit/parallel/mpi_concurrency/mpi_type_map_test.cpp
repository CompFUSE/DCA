// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner    (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests mpi_type_map.hpp.
// It is run with 1 MPI process.

#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

#include "gtest/gtest.h"
#include <complex>

enum Enum { A, B, C };

using Types = ::testing::Types<Enum, bool, int, unsigned long int, char, float, double,
                               std::complex<double>, std::complex<float>>;

template <typename T>
class MPITypeMapTest : public ::testing::Test {};

TYPED_TEST_CASE(MPITypeMapTest, Types);

TYPED_TEST(MPITypeMapTest, Types) {
  using ScalarType = TypeParam;

  EXPECT_EQ(sizeof(ScalarType), dca::parallel::MPITypeMap<ScalarType>::factor());
  EXPECT_EQ(MPI_CHAR, dca::parallel::MPITypeMap<ScalarType>::value());
}
