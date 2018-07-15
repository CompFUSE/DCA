// Copyright (C) 2012 Jakob Progsch, Václav Zeman
// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the methods of thead_pool.hpp

#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"

namespace dca {
namespace parallel {

ThreadPool::ThreadPool(size_t n_threads)
    : tasks_(n_threads), queue_mutex_(n_threads), condition_(n_threads), stop_(false), active_id_(0) {
  // Start a loop on a different thread.
  for (size_t id = 0; id < n_threads; ++id)
    workers_.emplace_back(&ThreadPool::workerLoop, this, id);
}

ThreadPool::ThreadPool() : ThreadPool(std::thread::hardware_concurrency()) {}

ThreadPool::~ThreadPool() {
  stop_ = true;
  for (auto& condition : condition_)
    condition.notify_one();
  for (std::thread& worker : workers_)
    worker.join();
}

void ThreadPool::workerLoop(int id) {
  while (true) {
    std::packaged_task<void()> task;

    {  // Acquire new task.
      std::unique_lock<std::mutex> lock(queue_mutex_[id]);
      condition_[id].wait(lock, [this, id] { return stop_ || !tasks_[id].empty(); });
      // If all the work is done and no more will be posted, return.
      if (stop_ && tasks_[id].empty())
        return;
      task = std::move(tasks_[id].front());
      tasks_[id].pop();
    }
    // Execute task.
    task();
  }
}

}  // parallel
}  // dca
