// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests elementary_vector_operations.hpp.

#include "math_library/geometry_library/vector_operations/elementary_vector_operations.hpp"
#include <vector>
#include "gtest/gtest.h"

TEST(elementary_vector_operations, IS_LARGER_VECTOR) {
  std::vector<double> v1 = {0., 1.};
  std::vector<double> v2 = {0., 2.};
  std::vector<double> v3 = {1., 0.};
  std::vector<double> v1_copy = v1;

  EXPECT_TRUE(VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v2));
  EXPECT_FALSE(VECTOR_OPERATIONS::IS_LARGER_VECTOR(v2, v1));
  EXPECT_TRUE(VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v3));
  EXPECT_FALSE(VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v1_copy));
}
