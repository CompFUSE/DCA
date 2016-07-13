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

#ifndef DCA_CONCURRENCY_PARALLELIZATION_PTHREADS_H
#define DCA_CONCURRENCY_PARALLELIZATION_PTHREADS_H

#include "dca/concurrency/parallelization_template.h"
#include <pthread.h>
#include "dca/concurrency/interfaces/processor_grouping_interface_pthreads.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class parallelization<POSIX_LIBRARY> {
public:
  parallelization();
  ~parallelization();

  void execute(int N, void* (*start_routine)(void*), void* arg);

  template <typename domain_t>
  static std::pair<int, int> get_bounds(int id, int N, domain_t& dmn);

  static std::pair<int, int> get_bounds(int id, int N, std::pair<int, int>& current_bounds);

private:
  void fork(int N, void* (*start_routine)(void*), void* arg);
  void join();

private:
  processor_grouping<POSIX_LIBRARY> group;
};

parallelization<POSIX_LIBRARY>::parallelization() {}

parallelization<POSIX_LIBRARY>::~parallelization() {}

void parallelization<POSIX_LIBRARY>::execute(int N, void* (*start_routine)(void*), void* arg) {
  group.fork(N, start_routine, arg);
  group.join();
}

void parallelization<POSIX_LIBRARY>::fork(int N, void* (*start_routine)(void*), void* arg) {
  group.fork(N, start_routine, arg);
}

void parallelization<POSIX_LIBRARY>::join() {
  group.join();
}

template <typename domain_t>
std::pair<int, int> parallelization<POSIX_LIBRARY>::get_bounds(int id, int N, domain_t& dmn) {
  std::pair<int, int> current_bounds(0, dmn.get_size());

  return get_bounds(id, N, current_bounds);
}

std::pair<int, int> parallelization<POSIX_LIBRARY>::get_bounds(int id, int N,
                                                               std::pair<int, int>& current_bounds) {
  long long size = current_bounds.second - current_bounds.first;

  long long bounds_first, bounds_second;

  long long CPU_id = static_cast<long long>(id);
  long long np = static_cast<long long>(N);

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
        (bounds.first >= 0 && bounds.second <= (current_bounds.second - current_bounds.first) &&
         bounds.first < bounds.second))) {
    std::cout << "error in " << __PRETTY_FUNCTION__ << "\n\n";
    std::cout << "CPU-id :: " << CPU_id << "\n";
    std::cout << "np     :: " << np << "\n";

    std::cout << "bounds.first  :: " << bounds.first << "\n";
    std::cout << "bounds.second :: " << bounds.second << "\n";

    throw std::logic_error(__FUNCTION__);
  }

  bounds.first += current_bounds.first;
  bounds.second += current_bounds.first;

  return bounds;
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_PARALLELIZATION_PTHREADS_H
