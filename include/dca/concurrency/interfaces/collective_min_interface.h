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

#ifndef DCA_CONCURRENCY_COLLECTIVE_MIN_INTERFACE_H
#define DCA_CONCURRENCY_COLLECTIVE_MIN_INTERFACE_H

#include "dca/concurrency/concurrency_types.hpp"
#include "dca/concurrency/interfaces/processor_grouping_interface.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class collective_min_interface {
public:
  collective_min_interface(processor_grouping<LIBRARY>& grouping_ref);
  ~collective_min_interface();

  template <typename scalar_type>
  void min(scalar_type& value);

private:
  processor_grouping<LIBRARY>& grouping;
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_min_interface<LIBRARY>::collective_min_interface(processor_grouping<LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_min_interface<LIBRARY>::~collective_min_interface() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_min_interface<LIBRARY>::min(scalar_type& /*value*/) {}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_COLLECTIVE_MIN_INTERFACE_H
