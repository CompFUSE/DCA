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

#ifndef DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_H
#define DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_H

#include <map>
#include <vector>
#include "dca/concurrency/concurrency_types.hpp"
#include "dca/concurrency/interfaces/processor_grouping_interface.h"
#include "comp_library/linalg/src/vector.h"
#include "comp_library/linalg/src/matrix.h"
#include "comp_library/function_library/function.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class collective_sum_interface {
public:
  collective_sum_interface(processor_grouping<LIBRARY>& grouping_ref);
  ~collective_sum_interface();

  template <typename scalar_type>
  void sum(scalar_type& value);

  template <typename scalar_type>
  void sum(std::vector<scalar_type>& m);

  template <typename scalartype>
  void sum(std::map<std::string, std::vector<scalartype>>& m);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f,
           FUNC_LIB::function<scalar_type, domain>& f_target);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f);

  template <typename some_type>
  void sum_and_average(some_type& obj, int size);

private:
  processor_grouping<LIBRARY>& grouping;
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_sum_interface<LIBRARY>::collective_sum_interface(processor_grouping<LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_sum_interface<LIBRARY>::~collective_sum_interface() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(scalar_type& /*value*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(std::vector<scalar_type>& /*m*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(std::map<std::string, std::vector<scalar_type>>& /*m*/) {
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& /*f*/,
                                            FUNC_LIB::function<scalar_type, domain>& /*f_target*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& /*f*/) {
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename some_type>
void collective_sum_interface<LIBRARY>::sum_and_average(some_type& /*obj*/, int /*size*/) {}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_H
