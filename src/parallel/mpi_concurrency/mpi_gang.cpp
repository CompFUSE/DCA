// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of mpi_gang.hpp.

#include "dca/parallel/mpi_concurrency/mpi_gang.hpp"

#include <cmath>

#include "dca/util/integer_division.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

MPIGang::MPIGang(const dca::parallel::MPIProcessorGrouping& group, int min_size) {
  min_size = std::min(min_size, group.get_size());


  const int n_gangs = group.get_size() / min_size;
  const int gang_id = std::min(group.get_id() / min_size, n_gangs - 1);

  MPI_Comm_split(group.get(), gang_id, 0, &communicator_);
  MPI_Comm_size(communicator_, &size_);
  MPI_Comm_rank(communicator_, &id_);
}

MPIGang::~MPIGang() {
  MPI_Comm_free(&communicator_);
}

}  // namespace parallel
}  // namespace dca
