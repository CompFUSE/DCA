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

template <bool complex_phase>
class MC_accumulator_data {
public:
  MC_accumulator_data() {
    initialize(0);
  }

  double& get_Gflop() {
    return gflop_;
  }

  auto get_accumulated_sign() const {
    return accumulated_sign_;
  }

  int get_number_of_measurements() const {
    return number_of_measurements_;
  }

  auto get_average_sign() const {
    return accumulated_sign_ / static_cast<double>(number_of_measurements_);
  }

  void initialize(int dca_iteration) {
    dca_iteration_ = dca_iteration;

    gflop_ = 0.;

    current_sign_ = complex_phase ? 0. : 1;
    accumulated_sign_ = 0;

    number_of_measurements_ = 0;
  }

protected:
  int dca_iteration_;

  double gflop_;

  using SignType = std::conditional_t<complex_phase, std::complex<double>, int>;
  SignType current_sign_;
  SignType accumulated_sign_;

  int number_of_measurements_;
};

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATOR_DATA_HPP
