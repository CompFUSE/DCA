// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements the thread task handler, that assigns a task (walker|accumulator|walker and
// accumulator) to each thread.
// For performance reasons walker and accumulator threads should be created alternately, i.e.
// walker, accumulator, w, a, w, a, ... . If the number of walkers and accumulators differ, the
// remaining threads are created at the end, e.g. for 4 walkers and 2 accumulators this means:
// w, a, w, a, w, w.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_THREAD_TASK_HANDLER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_THREAD_TASK_HANDLER_HPP

#include <cassert>
#include <string>
#include <vector>
#include <iostream>

namespace dca {
namespace phys {
namespace solver {

class ThreadTaskHandler {
public:
  ThreadTaskHandler(const int num_walkers, const int num_accumulators,
                    const bool shared_thread = false)
      : thread_tasks_(generateThreadTasksVec(num_walkers, num_accumulators, shared_thread)) {}

  // Prints all thread ids and the corresponding tasks (walker|accumulator|walker and accumulator).
  void print() const {
    for (int i = 0; i < thread_tasks_.size(); ++i)
      std::cout << "\t thread-id : " << i << "  -->   (" << thread_tasks_[i] << ")\n";
  }

  // Maps the walker id to the index of the walker's rng.
  int walkerIDToRngIndex(const int walker_id) const {
    assert(walker_id >= 0 && walker_id < thread_tasks_.size());
    assert(thread_tasks_[walker_id] == "walker");

    int rng_index = 0;

    for (int i = 0; i < walker_id; ++i) {
      if (thread_tasks_[i] == "walker")
        ++rng_index;
    }
    return rng_index;
  }

  // Maps the thread id to the id of the accumulator.
  int IDToAccumIndex(const int thread_id) const {
    assert(thread_id >= 0 && thread_id < thread_tasks_.size());
    assert(thread_tasks_[thread_id] == "accumulator");

    int rng_index = 0;

    for (int i = 0; i < thread_id; ++i) {
      if (thread_tasks_[i] == "accumulator")
        ++rng_index;
    }
    return rng_index;
  }

  std::size_t size() const {
    return thread_tasks_.size();
  }

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
                                                         bool shared_thread) {
    std::vector<std::string> thread_tasks;

    if (shared_thread && num_walkers != num_accumulators) {
      std::cerr << "Shared walker and accumulator thread requested, but the number of accumulators "
                   "and walkers differ.\nDiabling shared thread.\n";
      shared_thread = false;
    }

    if (!shared_thread) {
      thread_tasks = std::vector<std::string>(num_walkers + num_accumulators, "undefined");
      const int min_num_walkers_num_accumulators = std::min(num_walkers, num_accumulators);

      for (int i = 0; i < 2 * min_num_walkers_num_accumulators; i += 2) {
        thread_tasks[i] = "walker";
        thread_tasks[i + 1] = "accumulator";
      }

      for (int i = 2 * min_num_walkers_num_accumulators; i < num_walkers + num_accumulators; ++i) {
        if (min_num_walkers_num_accumulators != num_walkers)
          thread_tasks[i] = "walker";
        else
          thread_tasks[i] = "accumulator";
      }
    }
    else {
      thread_tasks = std::vector<std::string>(num_walkers, "walker and accumulator");
    }

    return thread_tasks;
  }

  const std::vector<std::string> thread_tasks_;
};

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_QMCI_THREAD_TASK_HANDLER_HPP
