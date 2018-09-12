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

constexpr char input_file[] = INPUT_DIR "input_3x3.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TpAccumulatorGpuTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

uint loop_counter = 0;

 TEST_F(TpAccumulatorGpuTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{27, 24};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorGpuTest::BDmn::dmn_size(),
                                        TpAccumulatorGpuTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n);

  for (const dca::phys::FourPointType type :
       {dca::phys::PARTICLE_HOLE_TRANSVERSE, dca::phys::PARTICLE_HOLE_MAGNETIC,
        dca::phys::PARTICLE_HOLE_CHARGE, dca::phys::PARTICLE_PARTICLE_UP_DOWN}) {
    parameters_.set_four_point_type(type);

    dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::CPU> accumulatorHost(
        data_->G0_k_w_cluster_excluded, parameters_);
    dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::linalg::GPU> accumulatorDevice(
        data_->G0_k_w_cluster_excluded, parameters_);
    const int sign = 1;

    accumulatorDevice.resetAccumulation(loop_counter);
    accumulatorDevice.accumulate(M, config, sign);
    accumulatorDevice.finalize();

    accumulatorHost.resetAccumulation(loop_counter);
    accumulatorHost.accumulate(M, config, sign);
    accumulatorHost.finalize();

    const auto diff = dca::func::util::difference(accumulatorHost.get_sign_times_G4(),
                                                  accumulatorDevice.get_sign_times_G4());
    EXPECT_GT(5e-7, diff.l_inf);
    ++loop_counter;
  }
}

TEST_F(TpAccumulatorGpuTest, SumToAndFinalize) {
  dca::linalg::util::initializeMagma();

  parameters_.set_four_point_type(dca::phys::PARTICLE_HOLE_TRANSVERSE);

  using Accumulator =
      dca::phys::solver::accumulator::TpAccumulator<G0Setup::Parameters, dca::linalg::GPU>;
  Accumulator accumulator_sum(data_->G0_k_w_cluster_excluded, parameters_, 0);
  Accumulator accumulator1(data_->G0_k_w_cluster_excluded, parameters_, 1);
  Accumulator accumulator2(data_->G0_k_w_cluster_excluded, parameters_, 2);
  Accumulator accumulator3(data_->G0_k_w_cluster_excluded, parameters_, 3);

  auto prepare_configuration = [&](auto& M, auto& configuration, const auto& n) {
    ConfigGenerator::prepareConfiguration(M, configuration, TpAccumulatorGpuTest::BDmn::dmn_size(),
                                          TpAccumulatorGpuTest::RDmn::dmn_size(),
                                          parameters_.get_beta(), n);
  };

  const std::array<int, 2> n{3, 4};
  const int sign = -1;
  Sample M1, M2;
  Configuration config1, config2;
  prepare_configuration(config1, M1, n);
  prepare_configuration(config2, M2, n);

  accumulator_sum.resetAccumulation(loop_counter++);

  accumulator1.accumulate(M1, config1, sign);
  accumulator2.accumulate(M2, config2, sign);
  accumulator1.sumTo(accumulator_sum);
  accumulator2.sumTo(accumulator_sum);
  accumulator_sum.finalize();

  // Reset the G4 on the GPU to zero.
  accumulator3.resetAccumulation(loop_counter++);
  accumulator3.accumulate(M1, config1, sign);
  accumulator3.accumulate(M2, config2, sign);
  accumulator3.finalize();

  const auto diff = dca::func::util::difference(accumulator3.get_sign_times_G4(),
                                                accumulator_sum.get_sign_times_G4());
  EXPECT_GT(5e-7, diff.l_inf);
}
