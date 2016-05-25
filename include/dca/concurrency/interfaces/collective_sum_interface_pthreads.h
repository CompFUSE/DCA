// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_PTHREADS_H
#define DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_PTHREADS_H

#include "dca/concurrency/interfaces/collective_sum_interface.h"
#include <pthread.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <>
class collective_sum_interface<POSIX_LIBRARY> {
public:
  template <typename scalar_type>
  static void sum(scalar_type& value, scalar_type& result, pthread_mutex_t& mutex);

  template <typename scalar_type, class domain>
  static void sum(FUNC_LIB::function<scalar_type, domain>& f,
                  FUNC_LIB::function<scalar_type, domain>& f_result, pthread_mutex_t& mutex);

private:
};

template <typename scalar_type>
void collective_sum_interface<POSIX_LIBRARY>::sum(scalar_type& value, scalar_type& result,
                                                  pthread_mutex_t& mutex) {
  pthread_mutex_lock(&mutex);

  value += result;

  pthread_mutex_unlock(&mutex);
}

template <typename scalar_type, class domain>
void collective_sum_interface<POSIX_LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& f,
                                                  FUNC_LIB::function<scalar_type, domain>& f_result,
                                                  pthread_mutex_t& mutex) {
  pthread_mutex_lock(&mutex);

  for (int l = 0; l < f.size(); l++)
    f_result(l) += f(l);

  pthread_mutex_unlock(&mutex);
}

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_COLLECTIVE_SUM_INTERFACE_PTHREADS_H
