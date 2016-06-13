// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Base class containing bare data for a MC accumulator.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_DATA_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_DATA_HPP

namespace DCA {
namespace QMCI {
// DCA::QMCI::

class MC_accumulator_data {
public:
  double& get_Gflop();

  double& get_sign();
  double& get_number_of_measurements();

  void initialize(int dca_iteration);

protected:
  int DCA_iteration;

  double GFLOP;

  double current_sign;
  double accumulated_sign;

  double number_of_measurements;
};

double& MC_accumulator_data::get_Gflop() {
  return GFLOP;
}

double& MC_accumulator_data::get_sign() {
  return accumulated_sign;
}

double& MC_accumulator_data::get_number_of_measurements() {
  return number_of_measurements;
}

void MC_accumulator_data::initialize(int dca_iteration) {
  DCA_iteration = dca_iteration;

  GFLOP = 0.;

  current_sign = 1;
  accumulated_sign = 0;

  number_of_measurements = 0;
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_DATA_HPP
