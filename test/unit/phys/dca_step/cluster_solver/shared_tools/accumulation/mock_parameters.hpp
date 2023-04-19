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
#include "dca/util/type_utils.hpp"

namespace dca {
namespace testing {

template <class BaseTestSetup, typename AccumType, bool singleGSampling>
struct MockParameters {
public:
  using profiler_type = dca::profiling::NullProfiler;
  using MC_measurement_scalar_type = AccumType;
  static constexpr bool complex_g0 = dca::util::IsComplex_t<AccumType>::value;
  using RClusterDmn = typename BaseTestSetup::RDmn;
  using Scalar = AccumType;
  using Real = dca::util::RealAlias<AccumType>;
  
  
  double get_beta() const {
    return beta_;
  }

  int stamping_period() const {
    if constexpr (singleGSampling)
      return 1;
    else
      return 0;
  }
  double beta_;
};

// template <class BaseTestSetup, typename AccumType>
// int MockParameters<BaseTestSetup, AccumType, false>::stamping_period() const {
//   return 0;
// }

// template <class BaseTestSetup, typename AccumType>
// int MockParameters<BaseTestSetup, AccumType, true>::stamping_period() const {
//   return 1;
// }

}  // namespace testing
}  // namespace dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_MOCK_PARAMETERS_HPP
