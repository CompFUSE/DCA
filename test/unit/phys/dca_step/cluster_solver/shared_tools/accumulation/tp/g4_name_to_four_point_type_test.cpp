// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests g4_name_to_four_point_type.hpp.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_name_to_four_point_type.hpp"
#include "gtest/gtest.h"

using namespace dca::phys;

TEST(G4NameToFourPointTypeTest, G4nameToFourPointType) {
  using solver::accumulator::G4nameToFourPointType;

  EXPECT_EQ(PARTICLE_HOLE_TRANSVERSE, G4nameToFourPointType("G4-particle-hole-transverse"));
  EXPECT_EQ(PARTICLE_HOLE_MAGNETIC, G4nameToFourPointType("G4-particle-hole-magnetic"));
  EXPECT_EQ(PARTICLE_HOLE_CHARGE, G4nameToFourPointType("G4-particle-hole-charge"));
  EXPECT_EQ(PARTICLE_PARTICLE_UP_DOWN, G4nameToFourPointType("G4-particle-particle-up-down"));

  EXPECT_THROW(G4nameToFourPointType("G4-anything"), std::invalid_argument);
}

TEST(G4NameToFourPointTypeTest, FourPointTypeToG4name) {
  using solver::accumulator::FourPointTypeToG4name;

  EXPECT_EQ("G4-particle-hole-transverse", FourPointTypeToG4name(PARTICLE_HOLE_TRANSVERSE));
  EXPECT_EQ("G4-particle-hole-magnetic", FourPointTypeToG4name(PARTICLE_HOLE_MAGNETIC));
  EXPECT_EQ("G4-particle-hole-charge", FourPointTypeToG4name(PARTICLE_HOLE_CHARGE));
  EXPECT_EQ("G4-particle-particle-up-down", FourPointTypeToG4name(PARTICLE_PARTICLE_UP_DOWN));

  EXPECT_THROW(FourPointTypeToG4name(NONE), std::invalid_argument);
}
