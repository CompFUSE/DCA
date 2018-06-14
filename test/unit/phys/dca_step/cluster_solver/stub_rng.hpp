// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
