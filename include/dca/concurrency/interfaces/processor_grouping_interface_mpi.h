// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef DCA_CONCURRENCY_PROCESSOR_GROUPING_INTERFACE_MPI_H
#define DCA_CONCURRENCY_PROCESSOR_GROUPING_INTERFACE_MPI_H

#include "dca/concurrency/interfaces/processor_grouping_interface.h"
#include <cassert>
#include <mpi.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class processor_grouping<MPI_LIBRARY> {
public:
  processor_grouping();
  ~processor_grouping();

  void set();
  MPI_Comm get();

  int get_id();

  int get_Nr_threads();

  int first();
  int last();

private:
  int id;
  int Nr_threads;

  MPI_Comm MPI_communication;
};

processor_grouping<MPI_LIBRARY>::processor_grouping() : id(-1), Nr_threads(-1) {}

processor_grouping<MPI_LIBRARY>::~processor_grouping() {}

void processor_grouping<MPI_LIBRARY>::set() {
  MPI_communication = MPI_COMM_WORLD;

  MPI_Comm_size(MPI_COMM_WORLD, &Nr_threads);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
}

MPI_Comm processor_grouping<MPI_LIBRARY>::get() {
  return MPI_communication;
}

int processor_grouping<MPI_LIBRARY>::get_id() {
  assert(id > -1);
  return id;
}

int processor_grouping<MPI_LIBRARY>::get_Nr_threads() {
  assert(Nr_threads > -1);
  return Nr_threads;
}

int processor_grouping<MPI_LIBRARY>::first() {
  return 0;
}

int processor_grouping<MPI_LIBRARY>::last() {
  return Nr_threads - 1;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_PROCESSOR_GROUPING_INTERFACE_MPI_H
