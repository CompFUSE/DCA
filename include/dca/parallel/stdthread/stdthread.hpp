// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class provides an interface for parallelizing with Pthreads.
//
// TODO: Finish sum methods.

#ifndef DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
#define DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP

#include <iostream>
#include <thread>
#include <vector>

#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"

namespace dca {
namespace parallel {

class stdthread {
public:
  stdthread() = default;

  // Execute the function f(id, num_threads, args...) as num_threads asynchronous tasks with id in
  // [0, num_threads - 1]. Then wait for the completion of the tasks.
  template <class F, class... Args>
  void execute(int num_threads, F&& f, Args&&... args) {
    std::vector<std::future<void>> futures;
    auto& pool = ThreadPool::get_instance();
    pool.enlarge(num_threads);

    // Fork.
    for (int id = 0; id < num_threads; ++id)
      futures.emplace_back(
          pool.enqueue(std::forward<F>(f), id, num_threads, std::forward<Args>(args)...));
    // Join.
    for (auto& future : futures)
      future.wait();
  }

  // Returns \sum_{id = 0}^{num_threads -1} f(id, num_threads, args...).
  // Precondition: the return type of f can be initialized with 0.
  template <class F, class... Args>
  auto sumReduction(int num_threads, F&& f, Args&&... args) {
    using ReturnType = typename std::result_of<F(int, int, Args...)>::type;

    std::vector<std::future<ReturnType>> futures;
    auto& pool = ThreadPool::get_instance();
    pool.enlarge(num_threads);

    // Fork.
    for (int id = 0; id < num_threads; ++id)
      futures.emplace_back(
          pool.enqueue(std::forward<F>(f), id, num_threads, std::forward<Args>(args)...));
    // Reduce.
    ReturnType result = 0;
    for (auto& future : futures)
      result += future.get();

    return result;
  }

  friend std::ostream& operator<<(std::ostream& some_ostream, const stdthread& this_concurrency);

private:
  static constexpr char parallel_type_str_[] = "stdthread";
};

}  //  parallel
}  //  dca

#endif  // DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
