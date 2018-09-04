// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements mpi_concurrency.hpp.

#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

#include <iostream>

namespace dca {
namespace parallel {

MPIConcurrency::MPIConcurrency(int argc, char** argv) : MPIInitializer(argc, argv) {
  if (!MPIProcessorGrouping::isValid()) {
    std::string error_string;
    if (MPIProcessorGrouping::get_size() == MPIProcessorGrouping::get_world_size())
      error_string = "No process could execute a CUDA kernel.";
    else
      error_string = "Process " + std::to_string(MPIProcessorGrouping::get_world_id()) +
                     "could not execute a CUDA kernel.";

    throw(std::logic_error(error_string));
  }
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
