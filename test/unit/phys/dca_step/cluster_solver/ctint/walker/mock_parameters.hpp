// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
// Giovanni Balduzzi(gbalduzz@ethz.phys.ch)
//
// This file provides the minimum parameters for initializing the time domain.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_MOCK_PARAMETERS_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_MOCK_PARAMETERS_HPP

namespace dca {
namespace ctint {
namespace testing {
// dca::testing::

class MockParameters {
private:
  double beta_;
  int time_intervals_;

public:
  MockParameters(double b, int i) : beta_(b), time_intervals_(i) {}

  double get_beta() const {
    return beta_;
  }

  int get_sp_time_intervals() const {
    return time_intervals_;
  }
};

}  // testing
}  // ctint
}  // dca

#endif  //  TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_MOCK_PARAMETERS_HPP
