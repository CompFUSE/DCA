// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the various types (channels) of vertex measurements.

#ifndef DCA_PHYS_VERTEX_MEASUREMENT_TYPE_HPP
#define DCA_PHYS_VERTEX_MEASUREMENT_TYPE_HPP

namespace dca {
namespace phys {
// dca::phys::

enum VertexMeasurementType {
  NONE,
  PARTICLE_HOLE_TRANSVERSE,
  PARTICLE_HOLE_MAGNETIC,
  PARTICLE_HOLE_CHARGE,
  PARTICLE_PARTICLE_SUPERCONDUCTING
};

}  // phys
}  // dca

#endif  // DCA_PHYS_VERTEX_MEASUREMENT_TYPE_HPP
