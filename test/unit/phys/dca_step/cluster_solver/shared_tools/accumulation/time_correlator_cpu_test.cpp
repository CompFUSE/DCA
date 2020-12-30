// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the TimeCorellator class.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "dca/linalg/util/handle_functions.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TimeCorrelatorTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_INT, input_file>;

using dca::linalg::CPU;

TEST_F(TimeCorrelatorTest, Accumulate) {
  dca::linalg::util::resizeHandleContainer(1);
  const int n_samples = 10;

  std::vector<Sample> M(n_samples);
  std::vector<Configuration> config(n_samples);

  for (int i = 0; i < n_samples; ++i) {
    const int n = 2 * n_samples - 1;
    ConfigGenerator::prepareConfiguration(config[i], M[i], TimeCorrelatorTest::BDmn::dmn_size(),
                                          TimeCorrelatorTest::RDmn::dmn_size(),
                                          parameters_.get_beta(), std::array<int, 2>{n, n});
  }
  const int n_correlations = 10;
  parameters_.set_time_correlation_window(n_correlations);

  using Correlator = dca::phys::solver::TimeCorrelator<Parameters, double, CPU>;
  Correlator correlator_(parameters_, 0);
  dca::phys::solver::G0Interpolation<CPU, double> g0(
      dca::phys::solver::ctint::details::shrinkG0(data_->G0_r_t));
  Correlator::setG0(g0);

  for (int i = 0; i < n_samples; ++i) {
    correlator_.compute_G_r_t(M[i], config[i], 1);
  }

  for (auto& c : correlator_.getCorrelators()) {
    EXPECT_GT(1., c.computeAutocorrelationTime());
  }
}
