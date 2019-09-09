// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the conversion of FourPointType to and from string.

#include "dca/phys/four_point_type.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
// dca::phys::

FourPointType stringToFourPointType(const std::string& name) {
  if (name == "PARTICLE_PARTICLE_UP_DOWN")
    return PARTICLE_PARTICLE_UP_DOWN;
  else if (name == "PARTICLE_HOLE_TRANSVERSE")
    return PARTICLE_HOLE_TRANSVERSE;
  else if (name == "PARTICLE_HOLE_MAGNETIC")
    return PARTICLE_HOLE_MAGNETIC;
  else if (name == "PARTICLE_HOLE_CHARGE")
    return PARTICLE_HOLE_CHARGE;
  else if (name == "PARTICLE_HOLE_LONGITUDINAL_UP_UP")
    return PARTICLE_HOLE_LONGITUDINAL_UP_UP;
  else if (name == "PARTICLE_HOLE_LONGITUDINAL_UP_DOWN")
    return PARTICLE_HOLE_LONGITUDINAL_UP_DOWN;
  else
    throw std::logic_error("Invalid four point mode: " + name);
}

std::string toString(const FourPointType type) {
  switch (type) {
    case PARTICLE_PARTICLE_UP_DOWN:
      return "PARTICLE_PARTICLE_UP_DOWN";
    case PARTICLE_HOLE_TRANSVERSE:
      return "PARTICLE_HOLE_TRANSVERSE";
    case PARTICLE_HOLE_MAGNETIC:
      return "PARTICLE_HOLE_MAGNETIC";
    case PARTICLE_HOLE_CHARGE:
      return "PARTICLE_HOLE_CHARGE";
    case PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      return "PARTICLE_HOLE_LONGITUDINAL_UP_UP";
    case PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      return "PARTICLE_HOLE_LONGITUDINAL_UP_DOWN";
    default:
      throw std::logic_error("Invalid four point mode.");
  }
}

}  // namespace phys
}  // namespace dca
