// Copyright 2015 ETH Zurich.
// Description.
//
// Author: Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich

#ifndef DCA_MPI_TEST_ENVIRONMENT_HPP
#define DCA_MPI_TEST_ENVIRONMENT_HPP

#include <string>

class dca_mpi_test_environment : public ::testing::Environment {
public:
  using concurrency_type =  COMP_LIB::parallelization<COMP_LIB::MPI_LIBRARY>;

  dca_mpi_test_environment(int argc,char *argv[], std::string file_name)
    : concurrency(argc, argv)
    , input_file(file_name)
  {}

  concurrency_type concurrency;
  std::string      input_file;
};

#endif  // DCA_MPI_TEST_ENVIRONMENT_HPP
