// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements thread_task_handler.hpp.

#include "dca/phys/dca_step/cluster_solver/posix_qmci/thread_task_handler.hpp"
#include <cassert>
#include <iostream>

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

ThreadTaskHandler::ThreadTaskHandler(const int num_walkers, const int num_accumulators)
    : thread_tasks_(generateThreadTasksVec(num_walkers, num_accumulators)) {}

void ThreadTaskHandler::print() const {
  for (int i = 0; i < thread_tasks_.size(); ++i) {
    if (thread_tasks_[i] == "walker")
      std::cout << "\t pthread-id : " << i << "  -->   (walker)\n";
    else
      std::cout << "\t pthread-id : " << i << "  -->   (accumulator)\n";
  }
}

int ThreadTaskHandler::walkerIDToRngIndex(const int walker_id) const {
  assert(walker_id >= 0 && walker_id < thread_tasks_.size());
  assert(thread_tasks_[walker_id] == "walker");

  int rng_index = 0;

  for (int i = 0; i < walker_id; ++i) {
    if (thread_tasks_[i] == "walker")
      ++rng_index;
  }
  return rng_index;
}

std::vector<std::string> ThreadTaskHandler::generateThreadTasksVec(const int num_walkers,
                                                                   const int num_accumulators) {
  std::vector<std::string> thread_tasks(num_walkers + num_accumulators, "undefined");

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
  return thread_tasks;
}

}  // solver
}  // phys
}  // dca
