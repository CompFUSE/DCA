// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements the methods in affinity.hpp.

// GNU extensions are required for linux-specific features for querying affinity
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "dca/parallel/stdthread/thread_pool/affinity.hpp"

#include <sched.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

// Includes for MacOS workaround
#ifdef __APPLE__
#include <pthread.h>
#include <mach/thread_act.h>
#include <sys/sysctl.h>
#include <bitset>
#include <limits>
#endif  // __APPLE__

namespace dca {
namespace parallel {
// dca::parallel::

// Workaround for MacOS
// Reference: https://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
#ifdef __APPLE__

#define SYSCTL_THREAD_COUNT "machdep.cpu.thread_count"

typedef struct cpu_set {
  using count_t = uint32_t;
  count_t count;
} cpu_set_t;

static constexpr int CPU_SETSIZE = std::numeric_limits<cpu_set_t::count_t>::digits;

static inline void CPU_ZERO(cpu_set_t* cpu_set) {
  cpu_set->count = 0;
}

static inline void CPU_SET(int num, cpu_set_t* cpu_set) {
  cpu_set->count |= (1 << num);
}

static inline int CPU_ISSET(int num, cpu_set_t* cpu_set) {
  return (cpu_set->count & (1 << num));
}

static inline int CPU_COUNT(cpu_set_t* cpu_set) {
  std::bitset<CPU_SETSIZE> bs(cpu_set->count);
  return bs.count();
}

int sched_getaffinity(pid_t /*pid*/, size_t /*cpu_size*/, cpu_set_t* cpu_set) {
  cpu_set_t::count_t core_count = 0;
  size_t len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_THREAD_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    std::cout << "Error while getting CPU thread count: " << ret << std::endl;
    return -1;
  }

  cpu_set->count = 0;
  for (cpu_set_t::count_t i = 0; i < core_count; ++i) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

// Binds the calling thread to the first available logical core.
int sched_setaffinity(pid_t /*pid*/, size_t cpu_size, cpu_set_t* cpu_set) {
  int core = 0;

  for (core = 0; core < 8 * cpu_size; ++core) {
    if (CPU_ISSET(core, cpu_set))
      break;
  }

  // std::cout << "Binding to logical core " << core << std::endl;
  thread_affinity_policy_data_t policy = {core};
  mach_port_t mach_thread = pthread_mach_thread_np(pthread_self());
  thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY, (thread_policy_t)&policy, 1);

  return 0;
}

#endif  // __APPLE__

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
