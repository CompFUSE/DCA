// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific
// publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

#ifndef DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_ACCUMULATOR_HPP
#define DCA_UNIT_PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_MOCK_ACCUMULATOR_HPP

#include "mock_parameters.hpp"
#include "mock_walker.hpp"

class MockAccumulator {
public:
  using my_parameters_type = Parameters;
  using my_MOMS_type = MOMS_t;

  MockAccumulator() {}
  MockAccumulator(Parameters &p_ref, MOMS_t &m_ref, int id0);
  void update_from(const MockWalker &wlk) {
    fake_configuration = wlk.get_configuration();
  }
  void measure();
  void sum_to(MockAccumulator& other){other.n_measurments+=n_measurments;}

  void initialize(int &iteration){};
  void finalize() {}
  double get_Gflop() { return 0; }
  double get_sign() { return 0; }
  int get_number_of_measurements() { return n_measurments; }

protected:
  MOMS_t MOMS;
  my_parameters_type parameters;
  std::vector<int> &get_configuration();

private:
  std::vector<int> fake_configuration;
  int id = -1;

public:
  int n_measurments = 0;
};

MockAccumulator::MockAccumulator(Parameters &p_ref, MOMS_t &m_ref, int id0)
        : parameters(p_ref), MOMS(m_ref), id(id0), fake_configuration(1, 0) {}

void MockAccumulator::measure() {
  if (fake_configuration[0] != 1)
    throw std::logic_error("configuration not loaded before measurment");
  fake_configuration[0] = 0;
  n_measurments++;
}


std::vector<int> &MockAccumulator::get_configuration() {
  return fake_configuration;
}

#endif
