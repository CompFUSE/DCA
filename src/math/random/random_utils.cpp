// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements random_utils.hpp.

#include "dca/math/random/random_utils.hpp"
#include <stdexcept>

namespace dca {
namespace math {
namespace random {
namespace detail {
// dca::math::random::detail::

int getGlobalId(const int local_id, const int proc_id, const int num_procs) {
  if (local_id < 0)
    throw std::logic_error("Local ID is invalid.");

  if (proc_id < 0 || proc_id >= num_procs)  // This implicitly checks for num_procs > 0.
    throw std::logic_error("Process ID is invalid.");

  return local_id * num_procs + proc_id;
}

uint64_t generateSeed(const int global_id, const uint64_t offset) {
  if (global_id < 0)
    throw std::logic_error("Global ID is invalid.");

  return hash(global_id + offset);
}

uint64_t hash(uint64_t key) {
  key ^= (key >> 33);
  key *= 0xbf58476d1ce4e5b9;
  key ^= (key >> 33);
  key *= 0x94d049bb133111eb;
  key ^= (key >> 33);

  return key;
}

}  // detail
}  // random
}  // math
}  // dca
