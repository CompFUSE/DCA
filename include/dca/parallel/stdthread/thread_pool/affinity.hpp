// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// These methods allow to retrieve and set the list of cores the current thread has affinity to.

#ifndef DCA_PARALLEL_STDTHREAD_THREAD_POOL_AFFINITY_HPP
#define DCA_PARALLEL_STDTHREAD_THREAD_POOL_AFFINITY_HPP

#include <vector>

namespace dca {
namespace parallel {
// dca::parallel::

// Returns a list of cores id for which the calling thread has affinity.
// If the macro __linux__ is not defined, returns an empty vector.
std::vector<int> get_affinity();

// Sets the affinity list of the current thread.
// If the macro __linux__ is not defined, performs a no-op.
void set_affinity(const std::vector<int>& cores);

// Returns the number of visible hardware cores.
// If the macro __linux__ is not defined, returns std::hardware_concurrency().
int get_core_count();

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_STDTHREAD_THREAD_POOL_AFFINITY_HPP
