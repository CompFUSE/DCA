// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface for parallelizing with Pthreads.
//
// TODO: Finish sum methods.

#ifndef DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
#define DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP

#include <iostream>
#include <thread>
#include <vector>

#include "dca/parallel/util/threading_data.hpp"

namespace dca {
namespace parallel {

class stdthread {
public:
  stdthread() : threads_(0), data_(0) {}

  void execute(int num_threads, void* (*start_routine)(void*), void* arg) {
    fork(num_threads, start_routine, arg);
    join();
  }

  friend std::ostream& operator<<(std::ostream& some_ostream, const stdthread& this_concurrency);

private:
  static constexpr char parallel_type_str_[] = "stdthread";
  void fork(int num_threads, void* (*start_routine)(void*), void* arg) {
    threads_.clear();
    data_.resize(num_threads);

    for (int id = 0; id < num_threads; id++) {
      data_[id].id = id;
      data_[id].num_threads = num_threads;
      data_[id].arg = arg;
      //
      threads_.push_back(std::thread(start_routine, (void*)(&data_[id])));
    }
  }

  void join() {
    for (int id = 0; id < threads_.size(); id++) {
      threads_[id].join();
    }
    threads_.clear();
  }

  std::vector<std::thread> threads_;
  std::vector<ThreadingData> data_;
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
