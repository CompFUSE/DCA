// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_PTHREADS_H
#define DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_PTHREADS_H

#include "dca/concurrency/interfaces/processor_grouping_interface.h"
#include <vector>
#include <pthread.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

struct posix_data {
public:
  posix_data() : id(-1), nr_threads(-1), args(NULL) {}

public:
  int id;
  int nr_threads;

  void* args;
};

template <>
class processor_grouping<POSIX_LIBRARY> {
public:
  processor_grouping();
  ~processor_grouping();

  void fork(int N, void* (*start_routine)(void*), void* arg);
  void join();

private:
  std::vector<pthread_t> pthread_vector;
  std::vector<posix_data> data_vector;
};

processor_grouping<POSIX_LIBRARY>::processor_grouping() : pthread_vector(0), data_vector(0) {}

processor_grouping<POSIX_LIBRARY>::~processor_grouping() {}

void processor_grouping<POSIX_LIBRARY>::fork(int N, void* (*routine)(void*), void* arg) {
  pthread_vector.resize(N);
  data_vector.resize(N);

  for (int l = 0; l < pthread_vector.size(); l++) {
    data_vector[l].id = l;
    data_vector[l].nr_threads = N;

    data_vector[l].args = arg;

    pthread_create(&pthread_vector[l], NULL, routine, (void*)(&data_vector[l]));
  }
}

void processor_grouping<POSIX_LIBRARY>::join() {
  for (int l = 0; l < pthread_vector.size(); l++)
    pthread_join(pthread_vector[l], NULL);

  pthread_vector.resize(0);
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_PROCESSOR_GROUPING_INTERFACE_PTHREADS_H
