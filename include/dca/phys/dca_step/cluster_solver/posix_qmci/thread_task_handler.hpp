// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements the thread task handler, that assigns a task (walker|accumulator) to each
// thread.
// For performance reasons walker and accumulator threads should be created alternately, i.e.
// walker, accumulator, w, a, w, a, ... . If the number of walkers and accumulators differ, the
// remaining threads are created at the end, e.g. for 4 walkers and 2 accumulators this means:
// w, a, w, a, w, w.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP

#include <cassert>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace posixqmci {
// dca::phys::solver::posixqmci::

class ThreadTaskHandler {
public:
  ThreadTaskHandler(const int num_walkers, const int num_accumulators, bool shared_walk_and_acc = false);

  // Prints all thread ids and the corresponding tasks (walker|accumulator).
  void print() const;

  // Maps the walker id to the index of the walker's rng.
  int walkerIDToRngIndex(const int walker_id) const;

  const std::string& getTask(const int id) const {
    assert(id < thread_tasks_.size());
    return thread_tasks_[id];
  }

  const std::vector<std::string>& getThreadTasksVec() const {
    return thread_tasks_;
  }

private:
  static std::vector<std::string> generateThreadTasksVec(const int num_walkers,
                                                         const int num_accumulators,
                                                         bool shared_walk);

  const std::vector<std::string> thread_tasks_;
};

}  // posixqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP
