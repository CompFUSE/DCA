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
  std::queue<std::packaged_task<void()>> tasks_;

  // synchronization
  std::mutex queue_mutex_;
  std::condition_variable condition_;
  bool stop_;
};

// the constructor just launches some amount of workers_
inline ThreadPool::ThreadPool(size_t n_threads) : stop_(false) {
  for (size_t i = 0; i < n_threads; ++i)
    workers_.emplace_back([this] {
      while (true) {
        std::packaged_task<void()> task;

        {  // Acquire new task.
          std::unique_lock<std::mutex> lock(this->queue_mutex_);
          this->condition_.wait(lock, [this] { return this->stop_ || !this->tasks_.empty(); });
          // If all the work is done and no more will be posted, return.
          if (this->stop_ && this->tasks_.empty())
            return;
          task = std::move(this->tasks_.front());
          this->tasks_.pop();
        }
        // Execute task.
        task();
      }
    });
}

// add new work item to the pool
template <class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<typename std::result_of<F(Args...)>::type> {
  using return_type = typename std::result_of<F(Args...)>::type;

  auto task =
      std::packaged_task<return_type()>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

  std::future<return_type> res = task.get_future();
  {
    std::unique_lock<std::mutex> lock(queue_mutex_);

    // don't allow enqueueing after stopping the pool
    if (stop_)
      throw std::runtime_error("enqueue on stopped ThreadPool");

    tasks_.emplace(std::move(task));
  }
  condition_.notify_one();
  return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool() {
  {
    std::unique_lock<std::mutex> lock(queue_mutex_);
    stop_ = true;
  }
  condition_.notify_all();
  for (std::thread& worker : workers_)
    worker.join();
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_STDTHREAD_THREAD_POOL_THREAD_POOL_HPP
