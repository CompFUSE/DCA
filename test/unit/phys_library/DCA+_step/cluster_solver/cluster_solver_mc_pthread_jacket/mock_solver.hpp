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
// Description: Mock MC solver for the posix test.

#ifndef DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_SOLVER_HPP
#define DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_SOLVER_HPP

#include "mock_parameters.hpp"
#include "mock_walker.hpp"
#include "mock_accumulator.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_pthread_jacket/posix_qmci_cluster_solver.h"
#include <iostream>
#include <stdexcept>

class MockSolver {
public:
  using walker_type = MockWalker;
  using accumulator_type = MockAccumulator;
  using this_MOMS_type = MOMS_t;
  using this_parameters_type = Parameters;

public:
  MockSolver(this_parameters_type& p_ref, this_MOMS_type& m_ref, bool standalone = true);
  void initialize(int i) {}
  double finalize(long int i) {
    return 0;
  }

protected:
  void update_shell(int i, int N, int N_k);
  void compute_error_bars(int tot_meas) {}
  void sum_measurements(int tot_meas) {}
  void symmetrize_measurements() {}

  double total_time = 0;
  int DCA_iteration = 0;
  this_parameters_type& parameters;
  this_MOMS_type MOMS;
  this_parameters_type::concurrency_type concurrency;
  MockAccumulator accumulator;  // pthread accumulator measurements get summed here
};

MockSolver::MockSolver(this_parameters_type& p_ref, this_MOMS_type& m_ref, bool standalone)
    : parameters(p_ref), MOMS(m_ref), concurrency(parameters.get_concurrency()) {
  if (standalone == true)
    throw std::logic_error("MockSolver is not called from pthread jacket");
}

// empty implementation
void MockSolver::update_shell(int i, int N, int N_k) {}

#endif
