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

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"

#include <array>
#include <limits>
#include <vector>
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/util/difference.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

using SpAccumulatorSingleMeasGTest = dca::testing::AccumulationTest<Scalar, 1, 3, 128, true>;

using MatrixPair = SpAccumulatorSingleMeasGTest::Sample;
using Configuration = SpAccumulatorSingleMeasGTest::Configuration;
using Parameters = SpAccumulatorSingleMeasGTest::Parameters;

TEST_F(SpAccumulatorSingleMeasGTest, Accumulate) {
  const std::array<int, 2> n{31, 28};
  MatrixPair M;
  Configuration config;
  prepareConfiguration(config, M, n, 0);

  using HostAccumulator = dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::CPU>;
  using MFunctionHost = typename HostAccumulator::MFunction;
  HostAccumulator accumulatorHost(
      parameters_);

  using DeviceAccumulator = dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::GPU>;
  DeviceAccumulator accumulatorDevice(
      parameters_);
  using MFunctionDev = typename DeviceAccumulator::MFunction;
  
  const int sign = 1;
  accumulatorDevice.resetAccumulation();
  accumulatorDevice.accumulate(M, config, sign);
  const MFunctionDev& mfunc1_dev(accumulatorDevice.get_single_measurement_sign_times_MFunction());
  
  accumulatorHost.resetAccumulation();
  accumulatorHost.accumulate(M, config, sign);
  const MFunctionHost& mfunc1_host(accumulatorHost.get_single_measurement_sign_times_MFunction());

  const auto diffMFunc1 = dca::func::util::difference(mfunc1_dev, mfunc1_host);
  EXPECT_GT(500 * std::numeric_limits<Parameters::MC_measurement_scalar_type>::epsilon(), diffMFunc1.l_inf);

  prepareConfiguration(config, M, n, 10000);
  accumulatorDevice.accumulate(M, config, sign);
  const MFunctionDev& mfunc2_dev(accumulatorDevice.get_single_measurement_sign_times_MFunction());
  accumulatorHost.accumulate(M, config, sign);
  const MFunctionHost& mfunc2_host(accumulatorHost.get_single_measurement_sign_times_MFunction());
  
  const auto diffMFunc2 = dca::func::util::difference(mfunc2_dev, mfunc2_host);
  EXPECT_GT(500 * std::numeric_limits<Parameters::MC_measurement_scalar_type>::epsilon(), diffMFunc2.l_inf);
  accumulatorDevice.finalize();
  accumulatorHost.finalize();

  const MFunctionHost& mfunc_accum_host(accumulatorHost.get_sign_times_M_r_w());
  const MFunctionDev& mfunc_accum_dev(accumulatorDevice.get_sign_times_M_r_w());
  
  const auto diffAccum = dca::func::util::difference(mfunc_accum_host,mfunc_accum_dev);

  EXPECT_GT(500 * std::numeric_limits<Parameters::MC_measurement_scalar_type>::epsilon(), diffAccum.l_inf);

  const auto diffMFunc2AccumDev = dca::func::util::difference(mfunc2_dev, mfunc_accum_dev);
  EXPECT_GT(diffMFunc2AccumDev.l1, 10000 * std::numeric_limits<Parameters::MC_measurement_scalar_type>::epsilon()) << diffMFunc2AccumDev.l_inf;
  const auto diffMFunc2AccumHost = dca::func::util::difference(mfunc2_host, mfunc_accum_host); 
  EXPECT_GT(diffMFunc2AccumHost.l1, 10000 * std::numeric_limits<Parameters::MC_measurement_scalar_type>::epsilon());
}

