// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
// Giovanni Balduzzi(gbalduzz@ethz.phys.ch)
//
// Provides a simple 1D domain whose size can change.

#ifndef TEST_UNIT_FUNCTION_VARIABLE_DOMAIN_HPP
#define TEST_UNIT_FUNCTION_VARIABLE_DOMAIN_HPP

namespace dca {
namespace func {
namespace testing {
// dca::testing::

class VariableDomain {
public:
  void static initialize(int n) {
    size_ = n;
  }
  // Used by dmn_variadic to compute the size.
  int static get_size() {
    return size_;
  }

  using element_type = void;

private:
  static int size_;
};

int VariableDomain::size_ = 0;

}  // testing
}  // func
}  // dca

#endif  //  TEST_UNIT_FUNCTION_VARIABLE_DOMAIN_HPP
