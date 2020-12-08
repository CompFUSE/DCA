// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Testing environment for DCA++ tests that use MPI.

#ifndef DCA_TESTING_DCA_MPI_TEST_ENVIRONMENT_HPP
#define DCA_TESTING_DCA_MPI_TEST_ENVIRONMENT_HPP

#include <string>
#include "gtest/gtest.h"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

namespace dca {
namespace testing {
// dca::testing::

struct DcaMpiTestEnvironment : public ::testing::Environment {
  using ConcurrencyType = dca::parallel::MPIConcurrency;

  DcaMpiTestEnvironment(ConcurrencyType& con, std::string file_name)
      : concurrency(con), input_file_name(file_name) {}

  ConcurrencyType& concurrency;
  std::string input_file_name;
};

}  // testing
}  // dca

#endif  // DCA_TESTING_DCA_MPI_TEST_ENVIRONMENT_HPP
