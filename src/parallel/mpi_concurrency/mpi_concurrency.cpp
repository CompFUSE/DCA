// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements mpi_concurrency.hpp.

#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

#include <iostream>

namespace dca {
namespace parallel {

MPIConcurrency::MPIConcurrency(int argc, char** argv)
    : MPIPacking(grouping_),
      MPICollectiveMax(grouping_),
      MPICollectiveMin(grouping_),
      MPICollectiveSum(grouping_) {
  // INTERNAL: Consider moving MPI_Init inside the MPIProcessorGrouping class.
  int provided = 0;
  constexpr int required = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, required, &provided);
  if (provided < required)
    throw(std::logic_error("MPI does not provide adequate thread support."));

  grouping_.reset(new MPIProcessorGrouping);

  if (!grouping_->isValid()) {  // Exit only from this process.
    std::cerr << "Process " << grouping_->get_world_id() << " is not valid.\nExiting" << std::endl;
    MPI_Finalize();
    exit(0);
  }
}

MPIConcurrency::~MPIConcurrency() {
  grouping_.release();
  MPI_Finalize();
}

constexpr char MPIConcurrency::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const MPIConcurrency& c) {
  o << '\n'
    << "concurrency type:" << c.parallel_type_str_ << '\n'
    << "number of processors:" << c.number_of_processors() << '\n'
    << "grouping first:" << c.first() << '\n'
    << "grouping last::" << c.last();
  return o;
}

}  // parallel
}  // dca
