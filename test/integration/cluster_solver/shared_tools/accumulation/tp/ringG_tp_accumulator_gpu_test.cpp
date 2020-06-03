// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Weile Wei (wwei9@lsu.edu)
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
#include "test/integration/cluster_solver/shared_tools/accumulation/tp/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4_multitransfer.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using DistributedTpAccumulatorGpuTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

uint loop_counter = 0;

TEST_F(DistributedTpAccumulatorGpuTest, Accumulate) {
    dca::linalg::util::initializeMagma();

    const std::array<int, 2> n{27, 24};
    Sample M;
    Configuration config;
    ConfigGenerator::prepareConfiguration(config, M, DistributedTpAccumulatorGpuTest::BDmn::dmn_size(),
          DistributedTpAccumulatorGpuTest::RDmn::dmn_size(),
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

    accumulatorDevice.resetAccumulation(loop_counter);
    accumulatorDevice.accumulate(M, config, sign);
    accumulatorDevice.finalize();

    accumulatorHost.resetAccumulation(loop_counter);
    accumulatorHost.accumulate(M, config, sign);
    accumulatorHost.finalize();

    ++loop_counter;

    auto& concurrency = parameters_.get_concurrency();
    for (int channel = 0; channel < accumulatorDevice.get_sign_times_G4().size(); ++channel) {
        auto G4_gpu = accumulatorDevice.get_sign_times_G4()[channel];
        auto G4_cpu = accumulatorHost.get_sign_times_G4()[channel];
        concurrency_.localSum(G4_gpu, concurrency.first());
        concurrency_.localSum(G4_cpu, concurrency.first());
        if (concurrency.get_id() == 0 && channel == 0){
            const auto diff = dca::func::util::difference(G4_cpu, G4_gpu);
            EXPECT_GT(5e-7, diff.l_inf);
            EXPECT_GT(5e-7, diff.l1);
            EXPECT_GT(5e-7, diff.l2);
        }
    }
}
