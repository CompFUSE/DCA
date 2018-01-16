// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides aa a testing util a mock rng with predetermined values.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_MOCK_RNG_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_MOCK_RNG_HPP

#include <vector>
#include <stdexcept>

namespace dca {
namespace testing {
// dca::testing::

class MockRng {
public:
  MockRng(int a = 0, int b = 0, int c = 0, int d = 0);  // For compatibility with solver.

  MockRng(const std::vector<double>& values) : val_(values) {}
  void setNewValues(const std::vector<double>& values) {
    val_ = values;
    index = 0;
  }
  double operator()() {
    if (index == val_.size()) {
      throw(std::out_of_range("no more values."));
    }
    return val_[index++];
  }

private:
  std::vector<double> val_;
  int index = 0;
};

MockRng::MockRng(int, int, int, int) {
  // Threat the rng parameters values as don't care.
  val_ = std::vector<double>(200,0.5);
}

}  // testing
}  // dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_MOCK_RNG_HPP
