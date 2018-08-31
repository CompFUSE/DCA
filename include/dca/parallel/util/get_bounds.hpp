// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides functions to distribute a domain among several threads.

#ifndef DCA_PARALLEL_UTIL_GET_BOUNDS_HPP
#define DCA_PARALLEL_UTIL_GET_BOUNDS_HPP

#include <utility>

namespace dca {
namespace parallel {
namespace util {
// dca::parallel::util::

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
std::pair<int, int> getBounds(int id, int num_threads, const Domain& dmn) {
  const std::pair<int, int> current_bounds(0, dmn.get_size());
  return getBounds(id, num_threads, current_bounds);
}

}  // util
}  // parallel
}  // dca

#endif  // DCA_PARALLEL_UTIL_GET_BOUNDS_HPP
