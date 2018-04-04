// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides as a testing util a stub rng with predetermined values.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_STUB_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_STUB_HPP

#include <vector>
#include <stdexcept>

namespace dca {
namespace testing {
// dca::testing::

class StubRng {
public:
  StubRng(int /*a*/ = 0, int /*b*/ = 0, int /*c*/ = 0, int /*d*/ = 0) {
    val_ = std::vector<double>(200, 0.5);
  }

  StubRng(const std::vector<double>& values) : val_(values) {}

  void setNewValues(const std::vector<double>& values) {
    val_ = values;
    index_ = 0;
  }

  double operator()() {
    if (index_ == val_.size()) {
      throw(std::out_of_range("No more values."));
    }
    return val_[index_++];
  }

private:
  std::vector<double> val_;
  int index_ = 0;
};

}  // testing
}  // dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_STUB_HPP
