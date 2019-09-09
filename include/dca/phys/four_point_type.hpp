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
// This file defines the types of four point functions.

#ifndef DCA_PHYS_FOUR_POINT_TYPE_HPP
#define DCA_PHYS_FOUR_POINT_TYPE_HPP

#include <string>

namespace dca {
namespace phys {
// dca::phys::

enum FourPointType : int {
  PARTICLE_HOLE_TRANSVERSE,
  PARTICLE_HOLE_MAGNETIC,
  PARTICLE_HOLE_CHARGE,
  PARTICLE_HOLE_LONGITUDINAL_UP_UP,
  PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
  PARTICLE_PARTICLE_UP_DOWN
};

FourPointType stringToFourPointType(const std::string& name);

std::string toString(const FourPointType type);

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_FOUR_POINT_TYPE_HPP
