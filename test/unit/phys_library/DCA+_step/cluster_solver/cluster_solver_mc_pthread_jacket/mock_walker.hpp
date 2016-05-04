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
// Description: Mock MC walker for the posix test.

#ifndef DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_WALKER_HPP
#define DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_WALKER_HPP

#include "mock_parameters.hpp"
using Parameters = MockParameters;
class MOMS_t {};

class MockWalker {
public:
  using rng_type = typename Parameters::random_number_generator;

public:
  MockWalker(Parameters& pars0, MOMS_t& moms, rng_type& rng0, int id0);
  void initialize() {}
  void do_sweep() {}
  void do_step() {}
  bool& is_thermalized() {
    return thermalized;
  }
  std::vector<int> get_configuration() const {
    return fake_configuration;
  }

  void static reset(int n_walk);
  static pthread_mutex_t walker_lock;
  static std::vector<bool> available_walker_ids;

private:
  bool thermalized = false;
  std::vector<int> fake_configuration;
  Parameters pars;
  rng_type rng;
  int thread_id;
};
pthread_mutex_t MockWalker::walker_lock;
std::vector<bool> MockWalker::available_walker_ids;

MockWalker::MockWalker(Parameters& pars0, MOMS_t& moms, rng_type& rng0, int id0)
    : pars(pars0), rng(rng0), thread_id(id0), fake_configuration(1, 1) {
  const int w_id = rng.id;
  std::cout << "Thread id: " << thread_id << "  Walker id:" << w_id << std::endl;
  const int nr_w = pars.get_nr_walkers();
  if (w_id >= nr_w or w_id < 0)
    throw std::logic_error("out of bond index");
  pthread_mutex_lock(&walker_lock);
  if (not available_walker_ids[w_id])
    throw std::logic_error("trying to use same id twice");
  available_walker_ids[w_id] = false;
  pthread_mutex_unlock(&walker_lock);
}

void MockWalker::reset(int n_w) {
  available_walker_ids = std::vector<bool>(n_w, true);
}

#endif
