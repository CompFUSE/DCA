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

#include "dca/config/threading.hpp"

namespace dca {
namespace parallel {

ThreadPool::ThreadPool(size_t n_threads) : stop_(false), active_id_(0) {
  core_count_ = get_core_count();
  master_affinity_ = get_affinity();

  enlarge(n_threads);
}

ThreadPool::~ThreadPool() {
  stop_ = true;
  for (auto& condition : condition_)
    condition->notify_one();
  for (std::thread& worker : workers_)
    worker.join();
}

void ThreadPool::enlarge(size_t n_threads) {
  if (workers_.size() >= n_threads)
    return;

  // Lock all queues to avoid modifications while enlarging the pool.
  std::vector<thread_traits::unique_lock> locks;
  for (auto& mutex_ptr : queue_mutex_)
    locks.emplace_back(thread_traits::unique_lock(*mutex_ptr));

  // Create enough resources for the new workers.
  queue_mutex_.reserve(n_threads);
  condition_.reserve(n_threads);
  tasks_.reserve(n_threads);

  // Release the locks on the existing queues.
  locks.clear();

  for (size_t id = workers_.size(); id < n_threads; ++id) {
    queue_mutex_.emplace_back(std::make_unique<thread_traits::mutex_type>());
    condition_.emplace_back(std::make_unique<thread_traits::condition_variable_type>());
    tasks_.emplace_back(std::make_unique<std::queue<std::packaged_task<void()>>>());

    // Start a loop on each new thread.
    workers_.emplace_back(&ThreadPool::workerLoop, this, id);
  }
}

void ThreadPool::workerLoop(int id) {
  // Set affinity.
  const int shift = core_count_ * id;
  std::vector<int> affinities;
  for (int x : master_affinity_)
    affinities.push_back(x + shift);

  set_affinity(affinities);

  while (true) {
    std::packaged_task<void()> task;

    {  // Acquire new task.
      thread_traits::unique_lock lock(*queue_mutex_[id]);
      condition_[id]->wait(lock, [this, id] { return stop_ || !tasks_[id]->empty(); });
      // If all the work is done and no more will be posted, return.
      if (stop_ && tasks_[id]->empty())
        return;
      task = std::move(tasks_[id]->front());
      tasks_[id]->pop();
    }
    // Execute task.
    task();
  }
}

}  // namespace parallel
}  // namespace dca
