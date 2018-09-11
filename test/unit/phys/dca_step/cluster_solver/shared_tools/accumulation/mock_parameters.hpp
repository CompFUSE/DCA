// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides a common setup for accumulation tests. It computes a mock configuration and
// M matrix over a single spin sector.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_MOCK_PARAMETERS_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_MOCK_PARAMETERS_HPP

#include "dca/profiling/null_profiler.hpp"

namespace dca {
namespace testing {

template <class BaseTestSetup, typename AccumType>
struct MockParameters {
public:
  using profiler_type = dca::profiling::NullProfiler;
  using MC_measurement_scalar_type = AccumType;

  using RClusterDmn = typename BaseTestSetup::RDmn;

  double get_beta() const {
    return beta_;
  }

  double beta_;
};

}  // testing
}  // dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_MOCK_PARAMETERS_HPP
