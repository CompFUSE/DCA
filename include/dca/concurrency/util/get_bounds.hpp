// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides functions to distribute a domain among several threads.

#ifndef DCA_PARALLEL_UTIL_GET_BOUNDS_HPP
#define DCA_PARALLEL_UTIL_GET_BOUNDS_HPP

#include <utility>

namespace dca {
namespace concurrency {
namespace util {
// dca::concurrency::util::

// Distributes the range defined by 'current_bounds' among 'num_threads' threads/processes and
// returns the bounds for 'id'.
//
// Preconditions: - id >= 0
//                - id < num_threads
//                - current_bounds.first <= current_bounds.second
std::pair<int, int> getBounds(int id, int num_threads, const std::pair<int, int>& current_bounds);

// Distributes the domain 'dmn' among 'num_threads' threads/processes and returns the bounds for
// 'id'.
//
// Preconditions: - id >= 0
//                - id < num_threads
//
// TODO: Add const to function parameter 'dmn'.
template <typename Domain>
std::pair<int, int> getBounds(int id, int num_threads, Domain& dmn) {
  const std::pair<int, int> current_bounds(0, dmn.get_size());
  return getBounds(id, num_threads, current_bounds);
}

}  // util
}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_UTIL_GET_BOUNDS_HPP
