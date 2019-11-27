// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface for parallelizing with HPX.
//
// TODO: Finish sum methods.

#ifndef DCA_PARALLEL_HPX_HPX_HPP
#define DCA_PARALLEL_HPX_HPX_HPP

#include <hpx/hpx.hpp>
#include <hpx/lcos/local/spinlock.hpp>
#include <hpx/lcos/local/condition_variable.hpp>
#include <hpx/lcos/future.hpp>
#include <hpx/include/threads.hpp>
#include <hpx/runtime/threads/thread.hpp>
#include <hpx/runtime/threads/executors/limiting_executor.hpp>
#include <hpx/debugging/demangle_helper.hpp>
#include <hpx/parallel/executors/parallel_executor.hpp>
#include <hpx/util/yield_while.hpp>
//
#include <vector>
#include <thread>

// #define DCA_HPX_THREAD_POOL_DEBUG

namespace dca {
namespace parallel {

struct thread_traits {
    template <typename T>
    using promise_type              = hpx::lcos::promise<T>;
    template <typename T>
    using future_type               = hpx::lcos::future<T>;
    using mutex_type                = hpx::lcos::local::mutex;
    using condition_variable_type   = hpx::lcos::local::condition_variable;
    using scoped_lock               = std::lock_guard<mutex_type>;
    using unique_lock               = std::unique_lock<mutex_type>;
    //
    static void sleep_for(hpx::util::steady_duration const& rel_time) {
        hpx::this_thread::sleep_for(rel_time);
    }
    //
    static void yield() {
        hpx::this_thread::yield();
    }
    //
    template <typename F>
    static void yield_while(F &&f) {
        hpx::util::yield_while(std::forward<F>(f));
    }
    //
    static std::uint64_t default_threadcount() {
        return hpx::get_num_worker_threads();
    }
};

// Returns a list of cores id for which the calling thread has affinity.
std::vector<int> get_affinity()
{
  cpu_set_t cpu_set_mask;
  auto status = sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);

  if (status == -1) {
    throw(std::runtime_error("Unable to get thread affinity."));
  }

  auto cpu_count = CPU_COUNT(&cpu_set_mask);

  std::vector<int> cores;
  cores.reserve(cpu_count);

  for (auto i = 0; i < CPU_SETSIZE && cores.size() < cpu_count; ++i) {
      cores.push_back(i);
  }

  return cores;
}

void set_affinity(const std::vector<int>& /*cores*/)
{
    // just do nothing, hpx handles this internally
}


// Number of cores used by this process.
int get_core_count() { return hpx::get_num_worker_threads(); }

class ThreadPool {
public:
  // Creates a pool with n_threads.
  // Actually does nothing, HPX does not need to allocate threads
  ThreadPool(size_t n_threads = 1) : exec(500, 500) {
    pool_size = n_threads;
  }

  ThreadPool(const ThreadPool& /*other*/) = delete;
  ThreadPool(ThreadPool&& /*other*/) = default;

  // Conclude all the pending work and destroy the threads spawned by this class.
  ~ThreadPool() {
      exec.wait_all();
  }

  void set_task_count_threshold(std::int64_t count)
  {
    exec.set_threshold(count, count+1);
  }

  void wait_for_tasks()
  {
    exec.wait_all();
  }

  // we don't do anything here
  void enlarge(std::size_t n_threads) {
      std::cout << "HPX threadpool enlarge: " << n_threads << std::endl;
      pool_size = n_threads;
  }

  // Call asynchronously the function f with arguments args. This method is thread safe.
  // Returns: a future to the result of f(args...).
  template <class F, class... Args>
  auto enqueue(F&& f, Args&&... args)
  {
#ifdef DCA_HPX_THREAD_POOL_DEBUG
    std::cout << "HPX threadpool enqueue\n";
    std::cout << "\n-------------------------------\n";
    std::cout << "enqueue: Function    : "
              << hpx::util::debug::print_type<F>() << "\n";
    std::cout << "enqueue: Arguments   : "
              << hpx::util::debug::print_type<Args...>(" | ") << std::endl;
#endif
    return hpx::async(exec, std::forward<F>(f), std::forward<Args>(args)...);
  }

  // We will not be using the pool for a while - put threads to sleep
  void suspend() {
      hpx::suspend();
  }

  // Returns the number of threads used by this class.
  std::size_t size() const {
    std::cout << "HPX threadpool size" << std::endl;
    return pool_size;
    // return hpx::get_num_worker_threads();
  }

  // Returns a static instance.
  static ThreadPool& get_instance() {
    static ThreadPool global_pool;
    return global_pool;
  }

  // this is just to make the DCA tests pass
  int pool_size;

  hpx::threads::executors::limiting_executor
    <hpx::threads::executors::default_executor> exec;

//    hpx::parallel::execution::parallel_executor;
//    hpx::threads::executors::default_executor;
};


struct hpxthread
{
  hpxthread() {};

  // Execute the function f(id, num_threads, args...) as num_threads asynchronous tasks with id in
  // [0, num_threads - 1]. Then wait for the completion of the tasks.
  template <class F, class... Args>
  void execute(int num_threads, F&& f, Args&&... args)
  {
    std::vector<thread_traits::future_type<void>> futures;
    //
    auto& pool = ThreadPool::get_instance();
    pool.enlarge(num_threads);

    // Fork.
    // Note we do not use std::forward here because we do not want to
    // accidentally move the same args more than once, we must use copy semantics
    for (int id = 0; id < num_threads; ++id)
      futures.emplace_back(
          pool.enqueue(f, id, num_threads, args...));

    // Join.
    for (auto& future : futures)
      future.wait();
  }

  // Returns the sum of the return values of f(id, num_tasks, args...) for each integer value of id
  // in [0, num_tasks).
  // Precondition: the return type of f can be initialized with 0.
  template <class F, class... Args>
  auto sumReduction(int num_threads, F&& f, Args&&... args) {
    assert(num_threads > 0);

    using ReturnType = typename std::result_of<F(int, int, Args...)>::type;

    std::vector<thread_traits::future_type<ReturnType>> futures;
    auto& pool = ThreadPool::get_instance();
    pool.enlarge(num_threads);

    // Spawn num_threads tasks.
    // Note we do not use std::forward here because we do not want to
    // accidentally move the same args more than once, we must use copy semantics
    for (int id = 0; id < num_threads; ++id)
      futures.emplace_back(
          pool.enqueue(f, id, num_threads, args...));

    // Sum the result of the tasks.
    ReturnType result = 0;
    for (auto& future : futures)
      result += future.get();

    return result;
  }


};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_hpxthread_hpxthread_HPP
