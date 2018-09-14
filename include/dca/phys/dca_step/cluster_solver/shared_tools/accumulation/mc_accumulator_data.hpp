// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Base class containing bare data for a MC accumulator.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATOR_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATOR_DATA_HPP

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

class MC_accumulator_data {
public:
  double& get_Gflop() {
    return GFLOP;
  }

  int get_accumulated_sign() const {
    return accumulated_sign;
  }

  double get_number_of_measurements() const {
    return number_of_measurements;
  }

  double get_average_sign() const {
    return static_cast<double>(accumulated_sign) / static_cast<double>(number_of_measurements);
  }

  void initialize(int dca_iteration) {
    DCA_iteration = dca_iteration;

    GFLOP = 0.;

    current_sign = 1;
    accumulated_sign = 0;

    number_of_measurements = 0;
  }

protected:
  int DCA_iteration;

  double GFLOP;

  int current_sign;
  int accumulated_sign;

  int number_of_measurements;
};

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATOR_DATA_HPP
