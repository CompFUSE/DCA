// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface to do a collective max with MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP

#include <memory>
#include <mpi.h>
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveMax {
public:
  MPICollectiveMax(const std::unique_ptr<const MPIProcessorGrouping>& grouping) : grouping_(grouping) {}

  template <typename Scalar>
  void max(Scalar& value) const {
    Scalar result;

    MPI_Allreduce(&value, &result, MPITypeMap<Scalar>::factor(), MPITypeMap<Scalar>::value(),
                  MPI_MAX, grouping_->get());

    value = result;
  }

private:
  const std::unique_ptr<const MPIProcessorGrouping>& grouping_;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MAX_HPP
