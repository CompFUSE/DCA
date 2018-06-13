// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
