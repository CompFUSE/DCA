// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the tetrahedron class.

#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"
#include "gtest/gtest.h"

TEST(TetrahedronTest, Constructor) {
  dca::math::geometry::tetrahedron<1> tetrahedron_1D;
  dca::math::geometry::tetrahedron<2> tetrahedron_2D;
  dca::math::geometry::tetrahedron<3> tetrahedron_3D;
}
