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

#ifndef PHYS_LIBRARY_VERTEX_MEASUREMENT_TYPE_HPP
#define PHYS_LIBRARY_VERTEX_MEASUREMENT_TYPE_HPP

enum VertexMeasurementType {
  NONE,
  PARTICLE_HOLE_TRANSVERSE,
  PARTICLE_HOLE_MAGNETIC,
  PARTICLE_HOLE_CHARGE,
  PARTICLE_PARTICLE_SUPERCONDUCTING
};

#endif  // PHYS_LIBRARY_VERTEX_MEASUREMENT_TYPE_HPP
