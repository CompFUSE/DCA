// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides methods to convert the names of G4 in specific channels to the corresponding
// four point types and vice versa.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_NAME_TO_FOUR_POINT_TYPE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_NAME_TO_FOUR_POINT_TYPE_HPP

#include <string>
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

FourPointType G4nameToFourPointType(const std::string& name) {
  if (name == "G4-particle-hole-transverse")
    return PARTICLE_HOLE_TRANSVERSE;
  if (name == "G4-particle-hole-magnetic")
    return PARTICLE_HOLE_MAGNETIC;
  if (name == "G4-particle-hole-charge")
    return PARTICLE_HOLE_CHARGE;
  if (name == "G4-particle-particle-up-down")
    return PARTICLE_PARTICLE_UP_DOWN;
  else
    throw std::invalid_argument("G4's name does not match any preset value");

  return NONE;
}

std::string FourPointTypeToG4name(const FourPointType channel) {
  if (channel == PARTICLE_HOLE_TRANSVERSE)
    return "G4-particle-hole-transverse";
  if (channel == PARTICLE_HOLE_MAGNETIC)
    return "G4-particle-hole-magnetic";
  if (channel == PARTICLE_HOLE_CHARGE)
    return "G4-particle-hole-charge";
  if (channel == PARTICLE_PARTICLE_UP_DOWN)
    return "G4-particle-particle-up-down";
  else
    throw std::invalid_argument("Four point type does not match any G4's name.");

  return "";
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_G4_NAME_TO_FOUR_POINT_TYPE_HPP
