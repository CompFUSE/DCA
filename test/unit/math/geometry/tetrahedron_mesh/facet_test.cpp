// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the facet class.

#include "dca/math/geometry/tetrahedron_mesh/facet.hpp"
#include "gtest/gtest.h"

TEST(FacetTest, Constructor) {
  dca::math::geometry::facet<2> facet_2D;
  dca::math::geometry::facet<3> facet_3D;
}
