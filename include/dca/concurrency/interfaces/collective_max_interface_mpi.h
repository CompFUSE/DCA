// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_CONCURRENCY_INTERFACES_COLLECTIVE_MAX_INTERFACE_MPI_H
#define DCA_CONCURRENCY_INTERFACES_COLLECTIVE_MAX_INTERFACE_MPI_H

#include "dca/concurrency/interfaces/collective_max_interface.h"
#include <mpi.h>
#include "dca/concurrency/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/concurrency/mpi_concurrency/mpi_type_map.hpp"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class collective_max_interface<MPI_LIBRARY> {
public:
  collective_max_interface(MPIProcessorGrouping& grouping_ref);
  ~collective_max_interface();

  template <typename scalar_type>
  void max(scalar_type& value);

private:
  MPIProcessorGrouping& grouping;
};

collective_max_interface<MPI_LIBRARY>::collective_max_interface(MPIProcessorGrouping& grouping_ref)
    : grouping(grouping_ref) {}

collective_max_interface<MPI_LIBRARY>::~collective_max_interface() {}

template <typename scalar_type>
void collective_max_interface<MPI_LIBRARY>::max(scalar_type& value) {
  scalar_type result;

  MPI_Allreduce(&value, &result, MPITypeMap<scalar_type>::factor(),
                MPITypeMap<scalar_type>::value(), MPI_MAX, grouping.get());

  value = result;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_COLLECTIVE_MAX_INTERFACE_MPI_H
