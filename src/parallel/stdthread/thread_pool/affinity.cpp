#include "dca/parallel/stdthread/thread_pool/affinity.hpp"

#include <iostream>
#include <cstdlib>

// GNU extensions are required for linux-specific features for querying affinity
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sched.h>

// returns a list of cores for which the calling thread has affinity
std::vector<int> get_affinity() {
  cpu_set_t cpu_set_mask;

  auto status = sched_getaffinity(0, sizeof(cpu_set_t), &cpu_set_mask);

  if (status == -1) {
    std::cerr << "error: unable to get affinity for thread" << std::endl;
    exit(1);
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
    std::cerr << "error: core count mismatch" << std::endl;
    exit(1);
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
