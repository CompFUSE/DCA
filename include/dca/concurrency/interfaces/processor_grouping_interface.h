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

#ifndef DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_H
#define DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_H

#include "dca/concurrency/concurrency_types.hpp"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class processor_grouping {
public:
  processor_grouping();
  ~processor_grouping();

  int get_id();

  int get_Nr_threads();

  int first();
  int last();

private:
  int id;          // This processors id within this grouping
  int Nr_threads;  // The number of processors within this grouping
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
processor_grouping<LIBRARY>::processor_grouping() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
processor_grouping<LIBRARY>::~processor_grouping() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
int processor_grouping<LIBRARY>::get_id() {
  return 0;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
int processor_grouping<LIBRARY>::get_Nr_threads() {
  return 1;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
int processor_grouping<LIBRARY>::first() {
  return 0;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
int processor_grouping<LIBRARY>::last() {
  return 0;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_H
