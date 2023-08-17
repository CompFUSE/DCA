// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file implements a comparison test for GPU vs host accumulation for a complex G0 case

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/platform/dca_gpu.h"

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator_gpu.hpp"

#include <array>
#include <limits>
#include <vector>

#include "dca/function/util/difference.hpp"

using Scalar = std::complex<double>;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

constexpr bool write_G4s = true;

#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
dca::parallel::MPIConcurrency* concurrency_ptr;
#else
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
dca::parallel::NoConcurrency* concurrency_ptr;
#endif

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_222-2_rashba.json";


using ConfigGenerator = dca::testing::AccumulationTest<std::complex<double>>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

template <typename SCALAR>
struct SpAccumulatorComplexG0GpuTest : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupBare<SCALAR, dca::testing::LatticeRashba,
                                            dca::ClusterSolverId::CT_AUX, input_file>;
  virtual void SetUp() {
    host_setup.SetUp();
    gpu_setup.SetUp();
  }

  virtual void TearDown() {}
  G0Setup host_setup;
  G0Setup gpu_setup;
};

//using Scalar = typename dca::config::McOptions::MCScalar;
using TestTypes = ::testing::Types<std::complex<double>>;
TYPED_TEST_CASE(SpAccumulatorComplexG0GpuTest, TestTypes);

#define TYPING_PREFACE                                            \
  using Scalar = TypeParam;                                       \
  using ConfigGenerator = dca::testing::AccumulationTest<Scalar>; \
  using Configuration = typename ConfigGenerator::Configuration;  \
  using Sample = typename ConfigGenerator::Sample;

TYPED_TEST(SpAccumulatorComplexG0GpuTest, Accumulate) {
  TYPING_PREFACE

  const std::array<int, 2> n{18, 22};

  using MatrixPair = Sample;
  //using Parameters = typename SpAccumulatorComplexG0GpuTest<Scalar>::Parameters;
  MatrixPair M;
  Configuration config;

  ConfigGenerator::prepareConfiguration(
      config, M, SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
      SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::RDmn::dmn_size(),
      this->host_setup.parameters_.get_beta(), n);

  dca::phys::solver::accumulator::SpAccumulator<decltype(this->host_setup.parameters_), dca::linalg::CPU> accumulatorHost(
      this->host_setup.parameters_);
  dca::phys::solver::accumulator::SpAccumulator<decltype(this->gpu_setup.parameters_), dca::linalg::GPU> accumulatorDevice(
      this->gpu_setup.parameters_);

  dca::util::SignType<Scalar> sign{1};
  accumulatorDevice.resetAccumulation();
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation();
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  const auto diff = dca::func::util::difference(accumulatorHost.get_sign_times_M_r_w(),
                                                accumulatorDevice.get_sign_times_M_r_w());

  EXPECT_GT(500 * std::numeric_limits<typename decltype(this->host_setup.parameters_)::Real>::epsilon(), diff.l_inf);
}

TYPED_TEST(SpAccumulatorComplexG0GpuTest, SumTo) {
  TYPING_PREFACE

  using Accumulator =
      dca::phys::solver::accumulator::SpAccumulator<decltype(this->gpu_setup.parameters_), dca::linalg::GPU>;

  Accumulator accumulator1(this->gpu_setup.parameters_);
  Accumulator accumulator2(this->gpu_setup.parameters_);
  Accumulator accumulator_sum(this->gpu_setup.parameters_);
  Accumulator accumulator3(this->gpu_setup.parameters_);

  const std::array<int, 2> n{3, 4};
  dca::util::SignType<Scalar> sign{-1};
  using MatrixPair = Sample;
  MatrixPair M1, M2;
  Configuration config1, config2;
  ConfigGenerator::prepareConfiguration(config1, M1, SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
					SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::RDmn::dmn_size(),  this->gpu_setup.parameters_.get_beta(), n);
  ConfigGenerator::prepareConfiguration(config2, M2, SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
					SpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::RDmn::dmn_size(),  this->gpu_setup.parameters_.get_beta(), n);

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

int main(int argc, char** argv) {
#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#else
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#endif

  dca::linalg::util::initializeMagma();

  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
