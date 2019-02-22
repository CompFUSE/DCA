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

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATION_TEST_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATION_TEST_HPP

#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/single_sector_accumulation_test.hpp"

#include <array>

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/mock_parameters.hpp"

namespace dca {
namespace testing {

template <typename AccumType, int n_bands = 2, int n_sites = 3, int n_frqs = 64>
class AccumulationTest : public SingleSectorAccumulationTest<AccumType, n_bands, n_sites, n_frqs> {
public:
  using TimeDmn = dca::func::dmn_0<dca::phys::domains::time_domain>;

  using BaseClass = SingleSectorAccumulationTest<AccumType, n_bands, n_sites, n_frqs>;
  using Parameters = MockParameters<BaseClass, AccumType>;
  using Configuration = std::array<typename BaseClass::Configuration, 2>;
  using Sample = std::array<typename BaseClass::Matrix, 2>;

protected:
  static void SetUpTestCase() {
    BaseClass::SetUpTestCase();

    // Initialize time domain.
    const int n_times = n_frqs;
    dca::phys::domains::time_domain::initialize(BaseClass::beta_, n_times);
  }

  void SetUp() {}

  void prepareConfiguration(Configuration& config, Sample& M, const std::array<int, 2> n) {
    for (int s = 0; s < 2; ++s)
      BaseClass::prepareConfiguration(config[s], M[s], n[s]);
  }

public:
  static void prepareConfiguration(Configuration& config, Sample& M, const int nb, const int nr,
                                   const double beta, const std::array<int, 2> n) {
    for (int s = 0; s < 2; ++s)
      BaseClass::prepareConfiguration(config[s], M[s], nb, nr, beta, n[s]);
  }

protected:
  Parameters parameters_{BaseClass::get_beta()};
};

}  // testing
}  // dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_ACCUMULATION_TEST_HPP
