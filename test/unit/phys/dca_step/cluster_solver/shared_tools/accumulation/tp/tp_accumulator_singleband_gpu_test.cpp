// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two particles accumulation on the GPU.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"

#include <array>
#include <functional>
#include <string>
#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TpAccumulatorGpuSinglebandTest =
    dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_AUX, input_file>;

unsigned int loop_id = 0;

TEST_F(TpAccumulatorGpuSinglebandTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{27, 24};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorGpuSinglebandTest::BDmn::dmn_size(),
                                        TpAccumulatorGpuSinglebandTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n);

  using namespace dca::phys;
  parameters_.set_four_point_channels(
      std::vector<FourPointType>{PARTICLE_HOLE_TRANSVERSE, PARTICLE_HOLE_MAGNETIC,
                                 PARTICLE_HOLE_CHARGE, PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                                 PARTICLE_HOLE_LONGITUDINAL_UP_DOWN, PARTICLE_PARTICLE_UP_DOWN});

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::CPU> accumulatorHost(
      data_->G0_k_w_cluster_excluded, parameters_);
  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::GPU> accumulatorDevice(
      data_->G0_k_w_cluster_excluded, parameters_);
  const int sign = 1;

  accumulatorDevice.resetAccumulation(0);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(0);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  for (std::size_t channel = 0; channel < accumulatorHost.get_sign_times_G4().size(); ++channel) {
    const auto diff = dca::func::util::difference(accumulatorHost.get_sign_times_G4()[channel],
                                                  accumulatorDevice.get_sign_times_G4()[channel]);
    EXPECT_GT(5e-7, diff.l_inf);
  }
}
