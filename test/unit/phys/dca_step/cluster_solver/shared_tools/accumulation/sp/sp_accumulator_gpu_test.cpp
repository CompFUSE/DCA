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

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"

#include <array>
#include "gtest/gtest.h"
#include <string>
#include <vector>

#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

using SpAccumulatorGpuTest = dca::testing::AccumulationTest<float, 1, 3, 128>;

using MatrixPair = SpAccumulatorGpuTest::Sample;
using Configuration = SpAccumulatorGpuTest::Configuration;
using Parameters = SpAccumulatorGpuTest::Parameters;

TEST_F(SpAccumulatorGpuTest, Accumulate) {
  const std::array<int, 2> n{31, 28};
  MatrixPair M;
  Configuration config;
  prepareConfiguration(config, M, n);

  dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::CPU> accumulatorHost(
      parameters_);
  dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::GPU> accumulatorDevice(
      parameters_);

  const int sign = 1;
  accumulatorDevice.resetAccumulation();
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation();
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  const auto diff = dca::func::util::difference(accumulatorHost.get_sign_times_M_r_w(),
                                                accumulatorDevice.get_sign_times_M_r_w());
  EXPECT_GT(5e-7, diff.l_inf);
}

TEST_F(SpAccumulatorGpuTest, SumTo) {
  using Accumulator = dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::GPU>;
  Accumulator accumulator1(parameters_);
  Accumulator accumulator2(parameters_);
  Accumulator accumulator_sum(parameters_);
  Accumulator accumulator3(parameters_);

  const std::array<int, 2> n{3, 4};
  const int sign = -1;
  MatrixPair M1, M2;
  Configuration config1, config2;
  prepareConfiguration(config1, M1, n);
  prepareConfiguration(config2, M2, n);

  accumulator1.accumulate(M1, config1, sign);
  accumulator2.accumulate(M2, config2, sign);
  accumulator1.sumTo(accumulator_sum);
  accumulator2.sumTo(accumulator_sum);
  accumulator_sum.finalize();

  accumulator3.accumulate(M1, config1, sign);
  accumulator3.accumulate(M2, config2, sign);
  accumulator3.finalize();

  const auto diff = dca::func::util::difference(accumulator3.get_sign_times_M_r_w(),
                                                accumulator_sum.get_sign_times_M_r_w());
  EXPECT_GT(5e-7, diff.l_inf);
}
