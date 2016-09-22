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

#ifndef DCA_CONCURRENCY_PARALLELIZATION_MPI_H
#define DCA_CONCURRENCY_PARALLELIZATION_MPI_H

#include "dca/concurrency/parallelization_template.h"
#include <mpi.h>
#include "dca/concurrency/mpi_concurrency/mpi_collective_max.hpp"
#include "dca/concurrency/mpi_concurrency/mpi_packing.hpp"
#include "dca/concurrency/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/concurrency/interfaces/collective_min_interface_mpi.h"
#include "dca/concurrency/interfaces/collective_sum_interface_mpi.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class parallelization<MPI_LIBRARY> : public MPIPacking,
                                     public MPICollectiveMax,
                                     public collective_min_interface<MPI_LIBRARY>,
                                     public collective_sum_interface<MPI_LIBRARY> {
public:
  parallelization(int argc, char* argv[]);
  ~parallelization();

  int id();

  int number_of_processors();

  int first();
  int last();

  template <typename object_type>
  bool broadcast(object_type& object, int root_id = 0);

  template <typename object_type>
  bool broadcast_object(object_type& object, int root_id = 0);

  template <typename domain_type>
  std::pair<int, int> get_bounds(domain_type& dmn);

private:
  MPIProcessorGrouping group;
};

parallelization<MPI_LIBRARY>::parallelization(int argc, char* argv[])
    : MPIPacking(group),
      MPICollectiveMax(group),
      collective_min_interface<MPI_LIBRARY>(group),
      collective_sum_interface<MPI_LIBRARY>(group) {
  MPI_Init(&argc, &argv);
  group.set();
}

parallelization<MPI_LIBRARY>::~parallelization() {
  MPI_Finalize();
}

int parallelization<MPI_LIBRARY>::first() {
  return group.first();
}

int parallelization<MPI_LIBRARY>::last() {
  return group.last();
}

int parallelization<MPI_LIBRARY>::id() {
  assert(group.get_id() > -1);
  return group.get_id();
}

int parallelization<MPI_LIBRARY>::number_of_processors() {
  assert(group.get_Nr_threads() > -1);
  return group.get_Nr_threads();
}

template <typename object_type>
bool parallelization<MPI_LIBRARY>::broadcast(object_type& object, int root_id) {
  assert(root_id > -1 and root_id < number_of_processors());

  int position = 0;

  if (id() == root_id) {
    int bufferSize = this->get_buffer_size(object);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, group.get());

    int* buffer = new int[bufferSize];

    this->pack(buffer, bufferSize, position, object);

    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, group.get());

    delete[] buffer;
  }
  else {
    int bufferSize(0);

    MPI_Bcast(&bufferSize, 1, MPI_INT, root_id, group.get());

    int* buffer = new int[bufferSize];

    // receive packed message
    MPI_Bcast(buffer, bufferSize, MPI_PACKED, root_id, group.get());

    this->unpack(buffer, bufferSize, position, object);

    delete[] buffer;
  }

  return true;
}

template <typename object_type>
bool parallelization<MPI_LIBRARY>::broadcast_object(object_type& object, int root_id) {
  assert(root_id > -1 and root_id < number_of_processors());

  int buffer_size = 0;

  if (id() == root_id) {
    buffer_size = object.get_buffer_size(*this);

    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, group.get());

    int* buffer = new int[buffer_size];

    int off_set = 0;
    object.pack(*this, buffer, buffer_size, off_set);

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, group.get());

    delete[] buffer;
  }
  else {
    MPI_Bcast(&buffer_size, 1, MPI_INT, root_id, group.get());

    int* buffer = new int[buffer_size];

    MPI_Bcast(buffer, buffer_size, MPI_PACKED, root_id, group.get());

    int off_set = 0;
    object.unpack(*this, buffer, buffer_size, off_set);

    delete[] buffer;
  }

  return true;
}

template <typename domain_type>
std::pair<int, int> parallelization<MPI_LIBRARY>::get_bounds(domain_type& dmn) {
  long long size = static_cast<long long>(dmn.get_size());

  long long bounds_first, bounds_second;

  long long CPU_id = static_cast<long long>(id());
  long long np = static_cast<long long>(number_of_processors());

  if (np < size) {
    bounds_first = (CPU_id * size) / np;
    bounds_second = ((CPU_id + 1) * size) / np;
  }
  else {
    if (CPU_id < size) {
      bounds_first = CPU_id;
      bounds_second = CPU_id + 1;
    }
    else {
      bounds_first = -1;
      bounds_second = -1;
    }
  }

  std::pair<int, int> bounds(static_cast<int>(bounds_first), static_cast<int>(bounds_second));

  if (!((bounds.first == -1 && bounds.second == -1) ||
        (bounds.first >= 0 && bounds.second <= dmn.get_size() && bounds.first < bounds.second))) {
    std::cout << "error in " << __PRETTY_FUNCTION__ << "\n\n";
    std::cout << "CPU-id :: " << CPU_id << "\n";
    std::cout << "np     :: " << np << "\n";

    std::cout << "bounds.first  :: " << bounds.first << "\n";
    std::cout << "bounds.second :: " << bounds.second << "\n";

    throw std::logic_error(__FUNCTION__);
  }

  return bounds;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_PARALLELIZATION_MPI_H
