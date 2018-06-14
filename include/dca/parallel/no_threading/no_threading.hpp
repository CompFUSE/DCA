// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class is the equivalent of Pthreading for serial execution.

#ifndef DCA_PARALLEL_NO_THREADING_HPP
#define DCA_PARALLEL_NO_THREADING_HPP

#include <iostream>
#include <stdexcept>
#include "dca/parallel/util/threading_data.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

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

  friend std::ostream& operator<<(std::ostream& some_ostream, const NoThreading& this_concurrency);

private:
  constexpr static char parallel_type_str_[] = "NoThreading";
  ThreadingData data_;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_THREADING_HPP
