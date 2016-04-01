// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests cluster_operations.hpp.

#include "cluster_operations.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "math_library/geometry_library/vector_operations/elementary_vector_operations.hpp"
#include <gtest/gtest.h>

template <typename T>
void printVectorSet(const std::vector<std::vector<T>>& vector_set,
                    const std::string name = "no-name") {
  std::cout << name << std::endl;
  for (const auto& vec : vector_set) {
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

class ClusterOperationsTest : public ::testing::Test {
protected:
  ClusterOperationsTest()
      : unsorted_set_{{0., 1.}, {-1., 1.}, {0., 0.}, {1., 0.}},
        sorted_set_(unsorted_set_),
        r_cluster_elements_{{-1.000000, 1.000000}, {-1.000000, 2.000000}, {0.000000, 0.000000},
                            {0.000000, 1.000000},  {0.000000, 2.000000},  {0.000000, 3.000000},
                            {1.000000, 1.000000},  {1.000000, 2.000000}},
        r_cluster_basis_{{1.00000, 0.00000}, {0.00000, 1.00000}},
        r_cluster_super_basis_{{2.00000, 2.00000}, {-2.00000, 2.00000}},
        k_cluster_elements_{{0.000000, 0.000000}, {0.000000, 3.141593}, {1.570796, 1.570796},
                            {1.570796, 4.712389}, {3.141593, 0.000000}, {3.141593, 3.141593},
                            {4.712389, 1.570796}, {4.712389, 4.712389}},
        k_cluster_basis_{{1.570796, 1.570796}, {-1.570796, 1.570796}},
        k_cluster_super_basis_{{6.283185, -0.000000}, {0.000000, 6.283185}} {
    // Vectors in the Brillouin zone are sorted according to VECTOR_OPERATIONS::IS_LARGER_VECTOR.
    std::sort(sorted_set_.begin(), sorted_set_.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR<double>);
  }
  std::vector<std::vector<double>> unsorted_set_;
  std::vector<std::vector<double>> sorted_set_;

  // 8-site cluster [2,]
  std::vector<std::vector<double>> r_cluster_elements_;
  std::vector<std::vector<double>> r_cluster_basis_;
  std::vector<std::vector<double>> r_cluster_super_basis_;

  std::vector<std::vector<double>> k_cluster_elements_;
  std::vector<std::vector<double>> k_cluster_basis_;
  std::vector<std::vector<double>> k_cluster_super_basis_;
};

TEST_F(ClusterOperationsTest, origin_index) {
  // Unsorted set (PARALLELEPIPEDUM)
  EXPECT_EQ(2, cluster_operations::origin_index(unsorted_set_, PARALLELLEPIPEDUM));

  // Sorted set (BRILLOUIN_ZONE)
  EXPECT_EQ(1, cluster_operations::origin_index(sorted_set_, BRILLOUIN_ZONE));
}

TEST_F(ClusterOperationsTest, index) {
  // Unsorted set (PARALLELEPIPEDUM)
  std::vector<double> vec1 = {-1., 1.};
  EXPECT_EQ(1, cluster_operations::index(vec1, unsorted_set_, PARALLELLEPIPEDUM));
  std::vector<double> vec2 = {-42., 24.};
  EXPECT_DEATH(cluster_operations::index(vec2, unsorted_set_, PARALLELLEPIPEDUM),
               "(index > -1 and index < elements.size())");

  // Sorted set (BRILLOUIN_ZONE)
  std::vector<double> vec3 = {0., 1.};
  EXPECT_EQ(2, cluster_operations::index(vec3, sorted_set_, BRILLOUIN_ZONE));
}

TEST_F(ClusterOperationsTest, translate_inside_cluster) {
  std::vector<double> input{1., 0.};
  std::vector<double> result{-1., 2.};
  EXPECT_EQ(result, cluster_operations::translate_inside_cluster(input, r_cluster_super_basis_));
}
