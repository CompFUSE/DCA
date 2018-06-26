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
// This file implements a thread pool based on the file
// https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h.
// See COPYING for the license of the original work.

#ifndef DCA_PARALLEL_STDTHREAD_THREAD_POOL_THREAD_POOL_HPP
#define DCA_PARALLEL_STDTHREAD_THREAD_POOL_THREAD_POOL_HPP

#include <condition_variable>
#include <iostream>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

namespace dca {
namespace parallel {

class ThreadPool {
public:
  ThreadPool(size_t n_threds);
  template <class F, class... Args>
  auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type>;
  ~ThreadPool();

  std::size_t size() const {
    return workers_.size();
  }

private:
  // need to keep track of threads so we can join them
  std::vector<std::thread> workers_;
  // the task queue
  std::vector<std::queue<std::packaged_task<void()>>> tasks_;

  // synchronization
  std::vector<std::mutex> queue_mutex_;
  std::vector<std::condition_variable> condition_;
  std::atomic<bool> stop_;
  unsigned int active_id_;
};

// the constructor just launches some amount of workers_
inline ThreadPool::ThreadPool(size_t n_threads)
    : tasks_(n_threads), queue_mutex_(n_threads), condition_(n_threads), stop_(false), active_id_(0) {
  for (size_t i = 0; i < n_threads; ++i)
    // Start a loop on a different thread.
    workers_.emplace_back(
        [this](const int id) {
          while (true) {
            std::packaged_task<void()> task;

            {  // Acquire new task.
              std::unique_lock<std::mutex> lock(this->queue_mutex_[id]);
              this->condition_[id].wait(
                  lock, [this, id] { return this->stop_ || !this->tasks_[id].empty(); });
              // If all the work is done and no more will be posted, return.
              if (this->stop_ && this->tasks_[id].empty())
                return;
              task = std::move(this->tasks_[id].front());
              this->tasks_[id].pop();
            }
            // Execute task.
            task();
          }
        },
        i);
}

// add new work item to the pool
template <class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<typename std::result_of<F(Args...)>::type> {
  using return_type = typename std::result_of<F(Args...)>::type;
  const int id = active_id_;

  auto task =
      std::packaged_task<return_type()>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task.get_future();
  {
    std::unique_lock<std::mutex> lock(queue_mutex_[id]);

    // don't allow enqueueing after stopping the pool
    if (stop_)
      throw std::runtime_error("enqueue on stopped ThreadPool");

    tasks_[id].emplace(std::move(task));
  }
  condition_[id].notify_one();

  active_id_ = (active_id_ + 1) % size();
  return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool() {
  stop_ = true;
  for (auto& condition : condition_)
    condition.notify_one();
  for (std::thread& worker : workers_)
    worker.join();
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_STDTHREAD_THREAD_POOL_THREAD_POOL_HPP
