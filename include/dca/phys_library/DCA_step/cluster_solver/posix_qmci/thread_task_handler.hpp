// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements the thread task handler, that assigns a task (walker|accumulator) to each
// thread.
// For performance reasons walker and accumulator threads should be created alternately, i.e.
// walker, accumulator, w, a, w, a, ... . If the number of walkers and accumulators differ, the
// remaining threads are created in the end, e.g. for 4 walkers and 2 accumulators this means: w, a,
// w, a, w, w.

#ifndef DCA_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP
#define DCA_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP

#include <string>
#include <vector>
#include <cassert>

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

class ThreadTaskHandler {
public:
  ThreadTaskHandler(const int num_walkers, const int num_accumulators);

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
                                                         const int num_accumulators);

  const int num_walkers_;
  const int num_accumulators_;
  const std::vector<std::string> thread_tasks_;
};

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_THREAD_TASK_HANDLER_HPP
