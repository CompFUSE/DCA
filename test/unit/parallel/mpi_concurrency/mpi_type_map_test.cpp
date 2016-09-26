// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests mpi_type_map.hpp.
// It is run with 1 MPI process.

#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"
#include "gtest/gtest.h"

using dca::parallel::MPITypeMap;

TEST(MPITypeMapTest, All) {
  EXPECT_EQ(1, MPITypeMap<bool>::factor());
  EXPECT_EQ(MPI_INT, MPITypeMap<bool>::value());

  EXPECT_EQ(1, MPITypeMap<char>::factor());
  EXPECT_EQ(MPI_CHAR, MPITypeMap<char>::value());

  EXPECT_EQ(1, MPITypeMap<int>::factor());
  EXPECT_EQ(MPI_INT, MPITypeMap<int>::value());

  EXPECT_EQ(1, MPITypeMap<std::size_t>::factor());
  EXPECT_EQ(MPI_UNSIGNED_LONG, MPITypeMap<std::size_t>::value());

  EXPECT_EQ(1, MPITypeMap<float>::factor());
  EXPECT_EQ(MPI_FLOAT, MPITypeMap<float>::value());

  EXPECT_EQ(1, MPITypeMap<double>::factor());
  EXPECT_EQ(MPI_DOUBLE, MPITypeMap<double>::value());

  EXPECT_EQ(2, MPITypeMap<std::complex<float>>::factor());
  EXPECT_EQ(MPI_FLOAT, MPITypeMap<std::complex<float>>::value());

  EXPECT_EQ(2, MPITypeMap<std::complex<double>>::factor());
  EXPECT_EQ(MPI_DOUBLE, MPITypeMap<std::complex<double>>::value());

  // Type without template specialization.
  struct MyStruct {
    int i;
    double d;
  };

  EXPECT_EQ(1, MPITypeMap<MyStruct>::factor());
  EXPECT_EQ(MPI_PACKED, MPITypeMap<MyStruct>::value());
}
