// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two particles accumulation on the GPU with
// the Rashba model.

#include "dca/config/profiler.hpp"
#include <complex>
using Scalar = std::complex<double>;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"

#include <array>
#include <functional>
#include <string>
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/util/difference.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_rashba.json";

using ConfigGenerator = dca::testing::AccumulationTest<std::complex<double>>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TpAccumulatorComplexG0GpuTest =
  dca::testing::G0Setup<Scalar, dca::testing::LatticeRashba, dca::ClusterSolverId::CT_AUX, input_file>;

uint loop_counter = 0;

TEST_F(TpAccumulatorComplexG0GpuTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{35, 0};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorComplexG0GpuTest::BDmn::dmn_size(),
                                        TpAccumulatorComplexG0GpuTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n);

  using namespace dca::phys;
  parameters_.set_four_point_channels(std::vector<FourPointType>{FourPointType::PARTICLE_PARTICLE_UP_DOWN});

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::CPU> accumulatorHost(
      data_->G0_k_w_cluster_excluded, parameters_);
  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::GPU> accumulatorDevice(
      data_->G0_k_w_cluster_excluded, parameters_);
  const int8_t sign = 1;

  accumulatorDevice.resetAccumulation(loop_counter);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(loop_counter);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  ++loop_counter;

  for (std::size_t channel = 0; channel < accumulatorHost.num_channels(); ++channel) {
    const auto diff = dca::func::util::difference(accumulatorHost.get_G4()[channel],
                                                  accumulatorDevice.get_G4()[channel]);
    EXPECT_GT(5e-7, diff.l_inf);
  }
}

