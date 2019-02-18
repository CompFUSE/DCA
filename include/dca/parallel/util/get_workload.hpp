// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a function that computes the local work size given local id and total amount
// of work.

#ifndef DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP
#define DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP

#include <cassert>

namespace dca {
namespace parallel {
namespace util {
// dca::parallel::util::

inline int getWorkload(const unsigned int total_work, const unsigned int n_workers,
                       const unsigned int id) {
  int work = total_work / n_workers;
  if (id < (total_work % n_workers))
    ++work;
  return work;
}

template <class Concurrency>
int getWorkload(const unsigned int total_work, const Concurrency& concurrency) {
  return getWorkload(total_work, concurrency.number_of_processors(), concurrency.id());
}

template <class Concurrency>
int getWorkload(const unsigned int total_work, const unsigned int n_local_workers,
                const unsigned int local_id, const Concurrency& concurrency) {
  assert(local_id < n_local_workers);
  const int local_work = getWorkload(total_work, concurrency);

  return getWorkload(local_work, n_local_workers, local_id);
}

}  // util
}  // parallel
}  // dca

#endif  // DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP
