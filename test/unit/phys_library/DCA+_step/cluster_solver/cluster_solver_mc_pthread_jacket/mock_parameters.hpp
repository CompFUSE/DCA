// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific
// publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Description: Mock set of parameters for the posix test.

#ifndef DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_PARAMETERS_HPP
#define DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_PARAMETERS_HPP
#undef MPI_SUPPORTED

#include <vector>
#include <pthread.h>
#include "comp_library/parallelization_library/include_parallelization_library.h"
#include "comp_library/profiler_library/profilers/null_profiler.h"
#include "mock_rng.hpp"

class MockParameters {
public:
  using profiler_type = PROFILER::no_profiling_mode;
  using concurrency_type = COMP_LIB::parallelization<COMP_LIB::SERIAL_LIBRARY>;
  using random_number_generator = MockRng;

public:
  MockParameters(int n_w = 2, int n_a = 2);
  ~MockParameters();
  int get_nr_walkers() {
    return nr_walkers;
  }
  int get_nr_accumulators() {
    return nr_accumulators;
  }
  int get_number_of_measurements() {
    return 10;
  }
  int get_warm_up_sweeps() {
    return 1;
  }
  int get_additional_steps() {
    return 1;
  }
  concurrency_type& get_concurrency();

private:
  int nr_accumulators;
  int nr_walkers;
  concurrency_type concurrency_obj;
  pthread_mutex_t parameters_lock;
};

MockParameters::MockParameters(int n_w, int n_a)
    : concurrency_obj(1, NULL), nr_walkers(n_w), nr_accumulators(n_a) {
  pthread_mutex_init(&parameters_lock, NULL);
}

MockParameters::~MockParameters() {
  pthread_mutex_destroy(&parameters_lock);
}

MockParameters::concurrency_type& MockParameters::get_concurrency() {
  return concurrency_obj;
}

#endif  // DCA_FAKE_PARAMETERS_HPP_H
