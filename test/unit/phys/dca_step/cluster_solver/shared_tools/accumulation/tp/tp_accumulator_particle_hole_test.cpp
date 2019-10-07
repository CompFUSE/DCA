// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the accumulation of the two-particle Green's function G4 in the particle-hole
// channels based on the relation between the particle-hole longitudinal up-up and up-down channels,
// and the particle-hole magnetic and charge channels.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <array>
#include <algorithm>  // for std::find_if
#include <limits>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_tp_accumulator_particle_hole_test.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;
using TpAccumulatorTest =
    dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_AUX, input_file>;

TEST_F(TpAccumulatorTest, ParticleHoleChannels) {
  Configuration config;
  Sample M;
  const std::array<int, 2> n{18, 22};
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorTest::BDmn::dmn_size(),
                                        TpAccumulatorTest::RDmn::dmn_size(), parameters_.get_beta(),
                                        n);

  parameters_.set_four_point_channels(std::vector<dca::phys::FourPointType>{
      dca::phys::PARTICLE_HOLE_LONGITUDINAL_UP_UP, dca::phys::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
      dca::phys::PARTICLE_HOLE_MAGNETIC, dca::phys::PARTICLE_HOLE_CHARGE});

  using TpAccumulatorType = dca::phys::solver::accumulator::TpAccumulator<Parameters>;
  TpAccumulatorType accumulator_ph(data_->G0_k_w_cluster_excluded, parameters_);

  const int sign = 1;
  accumulator_ph.accumulate(M, config, sign);

  accumulator_ph.finalize();

  const auto& G4 = accumulator_ph.get_sign_times_G4();

  // Check which element of the G4 container corresponds to which channel.
  const auto& G4_ph_magnetic =
      *std::find_if(G4.begin(), G4.end(), [](const TpAccumulatorType::TpGreensFunction& G4_channel) {
        return G4_channel.get_name() == "G4_PARTICLE_HOLE_MAGNETIC";
      });
  const auto& G4_ph_charge =
      *std::find_if(G4.begin(), G4.end(), [](const TpAccumulatorType::TpGreensFunction& G4_channel) {
        return G4_channel.get_name() == "G4_PARTICLE_HOLE_CHARGE";
      });
  const auto& G4_ph_long_up_up =
      *std::find_if(G4.begin(), G4.end(), [](const TpAccumulatorType::TpGreensFunction& G4_channel) {
        return G4_channel.get_name() == "G4_PARTICLE_HOLE_LONGITUDINAL_UP_UP";
      });
  const auto& G4_ph_long_up_down =
      *std::find_if(G4.begin(), G4.end(), [](const TpAccumulatorType::TpGreensFunction& G4_channel) {
        return G4_channel.get_name() == "G4_PARTICLE_HOLE_LONGITUDINAL_UP_DOWN";
      });

  // Compute G4 in the particle-hole longitudinal channels from the contributions in the p-h
  // magnetic and charge channels.
  TpAccumulatorType::TpGreensFunction G4_ph_long_up_up_check;
  TpAccumulatorType::TpGreensFunction G4_ph_long_up_down_check;

  for (int l = 0; l < G4_ph_long_up_up_check.size(); ++l) {
    G4_ph_long_up_up_check(l) = 0.5 * (G4_ph_charge(l) + G4_ph_magnetic(l));

    G4_ph_long_up_down_check(l) = 0.5 * (G4_ph_charge(l) - G4_ph_magnetic(l));
  }

  const auto diff_up_up = dca::func::util::difference(G4_ph_long_up_up, G4_ph_long_up_up_check);
  EXPECT_LT(diff_up_up.l_inf, 100 * std::numeric_limits<TpAccumulatorType::Real>::epsilon());

  const auto diff_up_down = dca::func::util::difference(G4_ph_long_up_down, G4_ph_long_up_down_check);
  EXPECT_LT(diff_up_down.l_inf, 100 * std::numeric_limits<TpAccumulatorType::Real>::epsilon());
}
