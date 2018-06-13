// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface to do a collective max with MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP

#include <mpi.h>
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveMax {
public:
  MPICollectiveMax(const MPIProcessorGrouping& grouping) : grouping_(grouping) {}

  template <typename Scalar>
  void max(Scalar& value) const {
    Scalar result;

    MPI_Allreduce(&value, &result, MPITypeMap<Scalar>::factor(), MPITypeMap<Scalar>::value(),
                  MPI_MAX, grouping_.get());

    value = result;
  }

private:
  const MPIProcessorGrouping& grouping_;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP
