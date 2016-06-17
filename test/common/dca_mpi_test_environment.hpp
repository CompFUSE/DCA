// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Testing environment for DCA++ tests that use MPI.

#ifndef DCA_TEST_COMMON_MPI_TEST_ENVIRONMENT_HPP
#define DCA_TEST_COMMON_MPI_TEST_ENVIRONMENT_HPP

#include <string>
#include "gtest/gtest.h"
#include "dca/concurrency/parallelization_mpi.h"

namespace dca {
namespace testing {
// dca::testing::

class DcaMpiTestEnvironment : public ::testing::Environment {
public:
  using ConcurrencyType = dca::concurrency::parallelization<dca::concurrency::MPI_LIBRARY>;

  DcaMpiTestEnvironment(int argc, char* argv[], std::string file_name)
      : concurrency(argc, argv), input_file_name(file_name) {}

  ConcurrencyType concurrency;
  std::string input_file_name;
};

}  // testing
}  // dca

#endif  // DCA_TEST_COMMON_MPI_TEST_ENVIRONMENT_HPP
