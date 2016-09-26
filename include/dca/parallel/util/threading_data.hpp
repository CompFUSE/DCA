// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This struct is used to create threads.

#ifndef DCA_PARALLEL_UTIL_THREADING_DATA_HPP
#define DCA_PARALLEL_UTIL_THREADING_DATA_HPP

#include <utility>

namespace dca {
namespace parallel {
// dca::parallel::

struct ThreadingData {
  ThreadingData() : id(-1), num_threads(-1), arg(nullptr) {}

  int id;
  int num_threads;
  void* arg;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_UTIL_THREADING_DATA
