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

#ifndef DCA_CONCURRENCY_INTERFACES_PRINT_ON_SHELL_INTERFACE_H
#define DCA_CONCURRENCY_INTERFACES_PRINT_ON_SHELL_INTERFACE_H

#include <iostream>
#include "dca/concurrency/concurrency_types.hpp"
#include "dca/concurrency/interfaces/processor_grouping_interface.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class print_on_shell_interface {
public:
  print_on_shell_interface(processor_grouping<LIBRARY>& grouping_ref);
  ~print_on_shell_interface();

  template <class whatever>
  void operator<<(whatever& whtvr);

  template <class whatever>
  void operator<<(const whatever& whtvr);

  template <class whatever>
  void operator<<(whatever* whtvr);

  template <class whatever>
  void operator<<(const whatever* whtvr);

private:
  processor_grouping<LIBRARY>& grouping;
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
print_on_shell_interface<LIBRARY>::print_on_shell_interface(processor_grouping<LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
print_on_shell_interface<LIBRARY>::~print_on_shell_interface() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <class whatever>
void print_on_shell_interface<LIBRARY>::operator<<(whatever& whtvr) {
  if (grouping.get_id() == grouping.first())
    std::cout << whtvr;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <class whatever>
void print_on_shell_interface<LIBRARY>::operator<<(const whatever& whtvr) {
  if (grouping.get_id() == grouping.first())
    std::cout << whtvr;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <class whatever>
void print_on_shell_interface<LIBRARY>::operator<<(whatever* whtvr) {
  if (grouping.get_id() == grouping.first())
    std::cout << whtvr;
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <class whatever>
void print_on_shell_interface<LIBRARY>::operator<<(const whatever* whtvr) {
  if (grouping.get_id() == grouping.first())
    std::cout << whtvr;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_PRINT_ON_SHELL_INTERFACE_H
