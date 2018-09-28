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
// This class provides a simple interface to the standard Message Passing Interface (MPI).

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <utility>
#include <stdexcept>

#include <mpi.h>

#include "dca/parallel/mpi_concurrency/mpi_collective_max.hpp"
#include "dca/parallel/mpi_concurrency/mpi_collective_min.hpp"
#include "dca/parallel/mpi_concurrency/mpi_collective_sum.hpp"
#include "dca/parallel/mpi_concurrency/mpi_initializer.hpp"
#include "dca/parallel/mpi_concurrency/mpi_packing.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/util/get_bounds.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPIConcurrency final : public virtual MPIInitializer,
                             private virtual MPIProcessorGrouping,
                             public MPIPacking,
                             public MPICollectiveMax,
                             public MPICollectiveMin,
                             public MPICollectiveSum {
public:
  MPIConcurrency(int argc, char** argv);

  inline int id() const {
    assert(MPIProcessorGrouping::get_id() > -1);
    return MPIProcessorGrouping::get_id();
  }
  int number_of_processors() const {
    assert(MPIProcessorGrouping::get_size() > -1);
    return MPIProcessorGrouping::get_size();
  }
  inline int first() const {
    return MPIProcessorGrouping::first();
  }
  inline int last() const {
    return MPIProcessorGrouping::last();
  }

  template <typename object_type>
  bool broadcast(object_type& object, int root_id = 0) const;

  template <typename object_type>
  bool broadcast_object(object_type& object, int root_id = 0) const;

  template <typename domain_type>
  std::pair<int, int> get_bounds(const domain_type& dmn) const {
    return util::getBounds(id(), number_of_processors(), dmn);
  }

  friend std::ostream& operator<<(std::ostream& some_ostream, const MPIConcurrency& this_concurrency);

private:
  constexpr static char parallel_type_str_[] = "MPIConcurrency";
};

template <typename object_type>
bool MPIConcurrency::broadcast(object_type& object, int root_id) const {
  assert(root_id > -1 and root_id < number_of_processors());

  int position = 0;

  if (id() == root_id) {
    int bufferSize = this->get_buffer_size(object);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, MPIProcessorGrouping::get());

    char* buffer = new char[bufferSize];

    this->pack(buffer, bufferSize, position, object);

    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, MPIProcessorGrouping::get());

    delete[] buffer;
  }
  else {
    int bufferSize(0);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, MPIProcessorGrouping::get());

    char* buffer = new char[bufferSize];

    // receive packed message
    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, MPIProcessorGrouping::get());

    this->unpack(buffer, bufferSize, position, object);

    delete[] buffer;
  }

  return true;
}

template <typename object_type>
bool MPIConcurrency::broadcast_object(object_type& object, int root_id) const {
  assert(root_id > -1 and root_id < number_of_processors());

  int buffer_size = 0;

  if (id() == root_id) {
    buffer_size = object.get_buffer_size(*this);

    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, MPIProcessorGrouping::get());

    char* buffer = new char[buffer_size];

    int off_set = 0;
    object.pack(*this, buffer, buffer_size, off_set);

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, MPIProcessorGrouping::get());

    delete[] buffer;
  }
  else {
    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, MPIProcessorGrouping::get());

    char* buffer = new char[buffer_size];

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, MPIProcessorGrouping::get());

    int off_set = 0;
    object.unpack(*this, buffer, buffer_size, off_set);

    delete[] buffer;
  }

  return true;
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_CONCURRENCY_HPP
