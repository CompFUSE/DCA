// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests four_point_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack, writing, and the computation of the 'exact'
//       momentum transfer and its index.

#include "dca/phys/parameters/four_point_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(FourPointParametersTest, DefaultValues) {
  dca::phys::params::FourPointParameters<2> pars;

  std::vector<double> momentum_transfer_input_check{0., 0.};

  EXPECT_EQ(dca::phys::NONE, pars.get_four_point_type());

  EXPECT_EQ(false, pars.accumulateG4ParticleHoleTransverse());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleMagnetic());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleCharge());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleLongitudinalUpUp());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleLongitudinalUpDown());
  EXPECT_EQ(false, pars.accumulateG4ParticleParticleUpDown());

  EXPECT_EQ(momentum_transfer_input_check, pars.get_four_point_momentum_transfer_input());
  EXPECT_EQ(0, pars.get_four_point_frequency_transfer());
  EXPECT_EQ(false, pars.compute_all_transfers());
}

TEST(FourPointParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::FourPointParameters<2> pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/four_point_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<double> momentum_transfer_input_check{3.14, -1.57};

  EXPECT_EQ(dca::phys::PARTICLE_PARTICLE_UP_DOWN, pars.get_four_point_type());

  EXPECT_EQ(true, pars.accumulateG4ParticleHoleTransverse());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleMagnetic());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleCharge());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleLongitudinalUpUp());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleLongitudinalUpDown());
  EXPECT_EQ(true, pars.accumulateG4ParticleParticleUpDown());

  EXPECT_EQ(momentum_transfer_input_check, pars.get_four_point_momentum_transfer_input());
  EXPECT_EQ(1, pars.get_four_point_frequency_transfer());
  EXPECT_EQ(true, pars.compute_all_transfers());
}

TEST(FourPointParametersTest, Setters) {
  dca::phys::params::FourPointParameters<2> pars;

  EXPECT_EQ(dca::phys::NONE, pars.get_four_point_type());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleTransverse());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleMagnetic());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleCharge());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleLongitudinalUpUp());
  EXPECT_EQ(false, pars.accumulateG4ParticleHoleLongitudinalUpDown());
  EXPECT_EQ(false, pars.accumulateG4ParticleParticleUpDown());

  pars.set_four_point_type(dca::phys::PARTICLE_HOLE_MAGNETIC);
  pars.accumulateG4ParticleHoleTransverse(true);
  pars.accumulateG4ParticleHoleMagnetic(true);
  pars.accumulateG4ParticleHoleCharge(true);
  pars.accumulateG4ParticleHoleLongitudinalUpUp(true);
  pars.accumulateG4ParticleHoleLongitudinalUpDown(true);
  pars.accumulateG4ParticleParticleUpDown(true);

  EXPECT_EQ(dca::phys::PARTICLE_HOLE_MAGNETIC, pars.get_four_point_type());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleTransverse());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleMagnetic());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleCharge());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleLongitudinalUpUp());
  EXPECT_EQ(true, pars.accumulateG4ParticleHoleLongitudinalUpDown());
  EXPECT_EQ(true, pars.accumulateG4ParticleParticleUpDown());
}

TEST(FourPointParametersTest, AccumulateG4) {
  dca::phys::params::FourPointParameters<2> pars;

  EXPECT_EQ(false, pars.accumulateG4());

  pars.set_four_point_type(dca::phys::PARTICLE_HOLE_MAGNETIC);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.set_four_point_type(dca::phys::NONE);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleHoleTransverse(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleHoleTransverse(false);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleHoleCharge(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleHoleCharge(false);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleHoleMagnetic(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleHoleMagnetic(false);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleHoleLongitudinalUpUp(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleHoleLongitudinalUpUp(false);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleHoleLongitudinalUpDown(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleHoleLongitudinalUpDown(false);
  EXPECT_EQ(false, pars.accumulateG4());

  pars.accumulateG4ParticleParticleUpDown(true);
  EXPECT_EQ(true, pars.accumulateG4());

  pars.accumulateG4ParticleParticleUpDown(false);
  EXPECT_EQ(false, pars.accumulateG4());
}
