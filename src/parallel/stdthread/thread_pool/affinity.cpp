// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// This file implements the methods in affinity.hpp.

#include "dca/parallel/stdthread/thread_pool/affinity.hpp"

#include <iostream>
#include <cstdlib>

// GNU extensions are required for linux-specific features for querying affinity
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sched.h>
#include <stdexcept>

namespace dca {
namespace parallel {
// dca::parallel::

std::vector<int> get_affinity() {
  cpu_set_t cpu_set_mask;

  auto status = sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);

  if (status == -1) {
    throw(std::runtime_error("Unable to get thread affinity."));
  }

  auto cpu_count = CPU_COUNT(&cpu_set_mask);

  std::vector<int> cores;
  cores.reserve(cpu_count);

  for (auto i = 0; i < CPU_SETSIZE && cores.size() < cpu_count; ++i) {
    if (CPU_ISSET(i, &cpu_set_mask)) {
      cores.push_back(i);
    }
  }

  if (cores.size() != cpu_count) {
    throw(std::logic_error("Core count mismatch."));
  }

  return cores;
}

void set_affinity(const std::vector<int>& cores) {
  cpu_set_t cpu_set_mask;
  CPU_ZERO(&cpu_set_mask);

  for (int core : cores) {
    CPU_SET(core, &cpu_set_mask);
  }

  sched_setaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);
}

int get_core_count() {
  cpu_set_t cpu_set_mask;
  sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);
  return CPU_COUNT(&cpu_set_mask);
}

}  // namespace parallel
}  // namespace dca
