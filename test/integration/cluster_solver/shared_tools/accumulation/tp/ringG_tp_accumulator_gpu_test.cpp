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

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_cpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"

#include <array>
#include <functional>
#include <string>
#include "gtest/gtest.h"
#include "dca/function/function.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/integration/cluster_solver/shared_tools/accumulation/tp/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR DCA_SOURCE_DIR "/test/integration/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4_multitransfer.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using DistributedTpAccumulatorGpuTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

uint loop_counter = 0;

TEST_F(DistributedTpAccumulatorGpuTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{24, 24};
  Sample M;
  Configuration config;
  unsigned long long unique_rng_skip = (n[0] * n[0] + 3 * n[0] + n[1] * n[1] + 3 * n[1]) * concurrency_.get_id() * 100;
  std::cout << "Rank " << concurrency_.get_id() << " is skipping " << unique_rng_skip << " random numbers\n";
  ConfigGenerator::prepareConfiguration(config, M, DistributedTpAccumulatorGpuTest::BDmn::dmn_size(),
                                        DistributedTpAccumulatorGpuTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n, unique_rng_skip);

  using namespace dca::phys;
  parameters_.set_four_point_channels(
      std::vector<FourPointType>{PARTICLE_HOLE_TRANSVERSE, PARTICLE_HOLE_MAGNETIC,
                                 PARTICLE_HOLE_CHARGE, PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                                 PARTICLE_HOLE_LONGITUDINAL_UP_DOWN, PARTICLE_PARTICLE_UP_DOWN});

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::GPU>
      accumulatorHost(data_->G0_k_w_cluster_excluded, parameters_);
  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::LINEAR, dca::linalg::GPU>
      accumulatorDevice(data_->G0_k_w_cluster_excluded, parameters_);
  const int sign = 1;

  accumulatorDevice.resetAccumulation(loop_counter);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(loop_counter);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  ++loop_counter;

  // auto& concurrency = parameters_.get_concurrency();

  if (concurrency_.get_id() == 0)
    std::cout << "\nCollecting Data from G4 distributed over " << concurrency_.number_of_processors()
              << "ranks\n";

  for (int channel = 0; channel < accumulatorDevice.get_G4().size(); ++channel) {
    DcaData<Parameters, dca::DistType::LINEAR>::TpGreensFunction& G4_gpu_dist =
        accumulatorDevice.get_nonconst_G4()[channel];
    DcaData<Parameters, dca::DistType::NONE>::TpGreensFunction& G4_gpu =
        accumulatorHost.get_nonconst_G4()[channel];
    concurrency_.sum(G4_gpu);
    G4_gpu_dist *= concurrency_.number_of_processors();
    const auto diff = dca::func::util::difference(G4_gpu, G4_gpu_dist);
    EXPECT_GT(5e-7, diff.l_inf);
    //EXPECT_GT(5e-7, diff.l1);
    //EXPECT_GT(5e-7, diff.l2);
  }
}
