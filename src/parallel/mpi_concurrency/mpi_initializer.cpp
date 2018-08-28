// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements mpi_initializer.hpp.

#include "dca/parallel/mpi_concurrency/mpi_initializer.hpp"

#include <stdexcept>
#include <mpi.h>

namespace dca {
namespace parallel {

MPIInitializer::MPIInitializer(int argc, char** argv) {
  int provided = 0;
  constexpr int required = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, required, &provided);
  if (provided < required)
    throw(std::logic_error("MPI does not provide adequate thread support."));
}

MPIInitializer::~MPIInitializer() {
  MPI_Finalize();
}

}  // parallel
}  // dca
