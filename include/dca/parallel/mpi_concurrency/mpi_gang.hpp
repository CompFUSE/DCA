// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class represent a subgroup of MPIProcessorGrouping.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_GANG_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_GANG_HPP

#include <mpi.h>

#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPIGang {
public:
  MPIGang(const MPIProcessorGrouping& group, int min_size);
  ~MPIGang();

  MPIGang(const MPIGang& other) = delete;

  int get_id() const {
    return id_;
  }
  int get_size() const {
    return size_;
  }

  auto get() const {
    return communicator_;
  }

private:
  MPI_Comm communicator_;
  int id_;
  int size_;
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_GANG_HPP
