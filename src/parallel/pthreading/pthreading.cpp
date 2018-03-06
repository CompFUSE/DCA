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
// This file implements pthreading.hpp

#include "dca/parallel/pthreading/pthreading.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

void Pthreading::execute(int num_threads, void* (*start_routine)(void*), void* arg) {
  fork(num_threads, start_routine, arg);
  join();
}

void Pthreading::fork(int num_threads, void* (*start_routine)(void*), void* arg) {
  pthreads_.resize(num_threads);
  data_.resize(num_threads);

  for (int id = 0; id < pthreads_.size(); id++) {
    data_[id].id = id;
    data_[id].num_threads = num_threads;
    data_[id].arg = arg;

    pthread_create(&pthreads_[id], NULL, start_routine, (void*)(&data_[id]));
  }
}

void Pthreading::join() {
  for (int id = 0; id < pthreads_.size(); id++)
    pthread_join(pthreads_[id], NULL);

  pthreads_.resize(0);
}

constexpr char Pthreading::concurrency_type_str_[];

std::ostream& operator<<(std::ostream& o, const Pthreading& c)
{
  o << '\n' << "concurrency type:" << c.concurrency_type_str_
    << '\n' << "number of posix threads:" << c.pthreads_.size();
  return o;
}

} // parallel
} // dca
