// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
