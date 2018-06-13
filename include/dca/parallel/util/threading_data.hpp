// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
