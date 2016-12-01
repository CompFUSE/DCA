// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests time_domain.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include <vector>
#include <gtest/gtest.h>

namespace dca {
namespace testing {
// dca::testing::

class MockParameters {
public:
  MockParameters(const double beta, const int sp_time_intervals)
      : beta_(beta), sp_time_intervals_(sp_time_intervals) {}

  double get_beta() const {
    return beta_;
  }
  int get_sp_time_intervals() const {
    return sp_time_intervals_;
  }

private:
  const double beta_;
  const int sp_time_intervals_;
};

}  // testing
}  // dca

TEST(TimeDomainTest, Initialize) {
  const double beta = 4.;
  const int sp_time_intervals = 2;
  const double eps = 1.e-10;

  dca::testing::MockParameters params(beta, sp_time_intervals);
  dca::phys::domains::time_domain::initialize(params);

  EXPECT_EQ(2 * (sp_time_intervals + 1), dca::phys::domains::time_domain::get_size());

  const std::vector<double> elements{-beta + eps, -beta / 2., -eps, eps, beta / 2., beta - eps};
  EXPECT_EQ(elements, dca::phys::domains::time_domain::get_elements());
}
