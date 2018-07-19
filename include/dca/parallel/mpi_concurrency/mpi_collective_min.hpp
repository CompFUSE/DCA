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
// This class provides an interface to do a collective min with MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MIN_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MIN_HPP

#include <memory>
#include <mpi.h>
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveMin {
public:
  MPICollectiveMin(const std::unique_ptr<const MPIProcessorGrouping>& grouping) : grouping_(grouping) {}

  template <typename Scalar>
  void min(Scalar& value) const {
    Scalar result;

    MPI_Allreduce(&value, &result, MPITypeMap<Scalar>::factor(), MPITypeMap<Scalar>::value(),
                  MPI_MIN, grouping_->get());

    value = result;
  }

private:
  const std::unique_ptr<const MPIProcessorGrouping>& grouping_;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_MIN_HPP
