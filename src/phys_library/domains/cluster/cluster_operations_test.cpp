// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests cluster_operations.h.

#include "cluster_operations.h"
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include "math_library/geometry_library/vector_operations/elementary_vector_operations.h"

template<typename T>
void print_elements(const std::vector<std::vector<T>>& elements, const std::string name="no-name") {
  std::cout << name << std::endl;
  for (const auto& vec : elements) {
    std::cout << "(";
    if (vec.size() > 0) {
      std::cout << vec[0];
      for (int i = 1; i < vec.size(); ++i) {
        std::cout << ", " << vec[i];
      }
    }
    std::cout << ")" << std::endl;
  }
}
  
TEST(cluster_operations, origin_index) {
  std::vector<std::vector<double>> unsorted_vec = {{0., 1.}, {-1., 1.}, {0.,0.}, {1., 0.}};
  std::vector<std::vector<double>> empty = {{}, {}};
  std::vector<std::vector<double>> one = {{1.,2.}};
  print_elements(unsorted_vec, "unsorted_vec");
  print_elements(empty, "empty");
  print_elements(one, "one");
  
  std::vector<std::vector<double>> sorted_vec(unsorted_vec);

  std::sort(sorted_vec.begin(), sorted_vec.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR<double>);
  print_elements(sorted_vec, "sorted_vec");
  
  std::vector<double> v1 = {0., 1.};
  std::vector<double> v2 = {0., 2.};
  std::vector<double> v3 = {1., 0.};
  std::cout << "VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v2) = " << VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v2) << std::endl;
  std::cout << "VECTOR_OPERATIONS::IS_LARGER_VECTOR(v2, v1) = " << VECTOR_OPERATIONS::IS_LARGER_VECTOR(v2, v1) << std::endl;
  std::cout << "VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v3) = " << VECTOR_OPERATIONS::IS_LARGER_VECTOR(v1, v3) << std::endl;

  EXPECT_EQ(1, cluster_operations::origin_index(sorted_vec, BRILLOUIN_ZONE));
  EXPECT_EQ(2, cluster_operations::origin_index(unsorted_vec, PARALLELLEPIPEDUM));
}
