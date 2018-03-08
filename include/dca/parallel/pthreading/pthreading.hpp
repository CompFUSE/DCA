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
// This class provides an interface for parallelizing with Pthreads.
//
// TODO: Finish sum methods.

#ifndef DCA_PARALLEL_PTHREADING_PTHREADING_HPP
#define DCA_PARALLEL_PTHREADING_PTHREADING_HPP

#include <vector>
#include <pthread.h>
#include <iostream>
#include "dca/parallel/util/threading_data.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class Pthreading {
public:
  Pthreading() : pthreads_(0), data_(0) {}

  void execute(int num_threads, void* (*start_routine)(void*), void* arg);

  // template <typename T>
  // static void sum(const T& value, T& result, pthread_mutex_t& mutex);
  // // TODO: Add const to function parameter 'f'.
  // template <typename T, typename Domain>
  // static void sum(func::function<T, Domain>& f, func::function<T, Domain>& f_result,
  //                 pthread_mutex_t& mutex);

  friend std::ostream& operator << (std::ostream& some_ostream, const Pthreading& this_concurrency);
private:
  constexpr static char parallel_type_str_[] = "PThreading";
  void fork(int num_threads, void* (*start_routine)(void*), void* arg);
  void join();

  std::vector<pthread_t> pthreads_;
  std::vector<ThreadingData> data_;
};

// template <typename T>
// void Pthreading::sum(const T& value, T& result, pthread_mutex_t& mutex) {
//   pthread_mutex_lock(&mutex);

//   // INTERNAL: Origininally it was "value += result;", but that must have been a typo.
//   result += value;

//   pthread_mutex_unlock(&mutex);
// }

// template <typename T, typename Domain>
// void Pthreading::sum(func::function<T, Domain>& f, func::function<T, Domain>& f_result,
//                      pthread_mutex_t& mutex) {
//   pthread_mutex_lock(&mutex);

//   for (int l = 0; l < f.size(); l++)
//     f_result(l) += f(l);

//   pthread_mutex_unlock(&mutex);
// }

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_PTHREADING_PTHREADING_HPP
