// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_CONCURRENCY_THREAD_MANAGER_SUM_H
#define DCA_CONCURRENCY_THREAD_MANAGER_SUM_H

#include <utility>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <typename concurrency_t>
class thread_manager_sum {
public:
  thread_manager_sum(concurrency_t& concurrency_ref);
  ~thread_manager_sum();

  template <typename domain_t>
  std::pair<int, int> get_bounds(domain_t& dmn);

  template <typename T>
  bool sum_and_check(T& result);

private:
  concurrency_t& concurrency_obj;
};

template <typename concurrency_t>
thread_manager_sum<concurrency_t>::thread_manager_sum(concurrency_t& concurrency_ref)
    : concurrency_obj(concurrency_ref) {}

template <typename concurrency_t>
thread_manager_sum<concurrency_t>::~thread_manager_sum() {}

template <typename concurrency_t>
template <typename domain_t>
std::pair<int, int> thread_manager_sum<concurrency_t>::get_bounds(domain_t& dmn) {
  return concurrency_obj.get_bounds(dmn);
}

template <typename concurrency_t>
template <typename T>
bool thread_manager_sum<concurrency_t>::sum_and_check(T& result) {
  concurrency_obj.sum(result);
  return true;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_THREAD_MANAGER_SUM_H
