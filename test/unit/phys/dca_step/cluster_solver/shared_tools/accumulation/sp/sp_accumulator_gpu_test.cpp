// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This file implements a comparison test for GPU vs host accumulation.

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/platform/dca_gpu.h"

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"

#include <array>
#include <limits>
#include <vector>

#include "dca/function/util/difference.hpp"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

//using Scalar = typename dca::config::McOptions::MCScalar;
template <typename SCALAR>
using SpAccumulatorGpuTest = dca::testing::AccumulationTest<SCALAR, 1, 3, 128>;
using TestTypes = ::testing::Types<float, double, std::complex<double>>;
TYPED_TEST_CASE(SpAccumulatorGpuTest, TestTypes);

TYPED_TEST(SpAccumulatorGpuTest, Accumulate) {
  using Scalar = TypeParam;
  using MatrixPair = typename SpAccumulatorGpuTest<Scalar>::Sample;
  using Configuration = typename SpAccumulatorGpuTest<Scalar>::Configuration;
  using Parameters = typename SpAccumulatorGpuTest<Scalar>::Parameters;
  const std::array<int, 2> n{31, 28};
  MatrixPair M;
  Configuration config;
  this->prepareConfiguration(config, M, n);

  dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::CPU> accumulatorHost(
      this->parameters_);
  dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::GPU> accumulatorDevice(
      this->parameters_);

  dca::util::SignType<Scalar> sign{1};
  accumulatorDevice.resetAccumulation();
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation();
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  const auto diff = dca::func::util::difference(accumulatorHost.get_sign_times_M_r_w(),
                                                accumulatorDevice.get_sign_times_M_r_w());

  EXPECT_GT(500 * std::numeric_limits<typename Parameters::Real>::epsilon(), diff.l_inf);
}

TYPED_TEST(SpAccumulatorGpuTest, SumTo) {
  using Scalar = TypeParam;
  using Parameters = typename SpAccumulatorGpuTest<Scalar>::Parameters;
  using Accumulator =
      dca::phys::solver::accumulator::SpAccumulator<Parameters, dca::linalg::GPU>;
  using Configuration = typename SpAccumulatorGpuTest<Scalar>::Configuration;

  Accumulator accumulator1(this->parameters_);
  Accumulator accumulator2(this->parameters_);
  Accumulator accumulator_sum(this->parameters_);
  Accumulator accumulator3(this->parameters_);

  const std::array<int, 2> n{3, 4};
  const int sign = -1;
  using MatrixPair = typename SpAccumulatorGpuTest<Scalar>::Sample;
  MatrixPair M1, M2;
  Configuration config1, config2;
  this->prepareConfiguration(config1, M1, n);
  this->prepareConfiguration(config2, M2, n);

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
