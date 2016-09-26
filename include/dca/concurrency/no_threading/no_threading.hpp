// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class is the equivalent of Pthreading for serial execution.

#ifndef DCA_PARALLEL_NO_THREADING_HPP
#define DCA_PARALLEL_NO_THREADING_HPP

#include "dca/concurrency/util/threading_data.hpp"

namespace dca {
namespace concurrency {
// dca::concurrency::

class NoThreading {
public:
  void execute(int num_threads, void* (*start_routine)(void*), void* arg) {
    for (int id = 0; id < num_threads; id++) {
      data_.id = id;
      data_.num_threads = num_threads;
      data_.arg = arg;
      start_routine(static_cast<void*>(&data_));
    }
  }

private:
  ThreadingData data_;
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_NO_THREADING_HPP
