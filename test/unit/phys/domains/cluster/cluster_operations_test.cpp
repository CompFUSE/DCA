// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests cluster_operations.hpp.

#include "dca/phys/domains/cluster/cluster_operations.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace testing {
// dca::testing::

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

}  // testing
}  // dca

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
    std::sort(sorted_set_.begin(), sorted_set_.end(), dca::math::util::isLessVector<double>);
  }
  std::vector<std::vector<double>> unsorted_set_;  // = {{0., 1.}, {-1., 1.}, {0., 0.}, {1., 0.}}
  // Vectors in the Brillouin zone are sorted according to dca::math::util::isLessVector.
  std::vector<std::vector<double>> sorted_set_;  // = {{-1., 1.}, {0., 0.}, {0., 1.}, {1., 0.}}

  // 8-site cluster {2,2}, {-2, 2}
  // Real space
  std::vector<std::vector<double>> r_cluster_elements_;
  std::vector<std::vector<double>> r_cluster_basis_;
  std::vector<std::vector<double>> r_cluster_super_basis_;
  // Momentum space
  std::vector<std::vector<double>> k_cluster_elements_;
  std::vector<std::vector<double>> k_cluster_basis_;
  std::vector<std::vector<double>> k_cluster_super_basis_;
};

using ClusterOperationsDeathTest = ClusterOperationsTest;

TEST_F(ClusterOperationsTest, origin_index) {
  // Unsorted set (PARALLELEPIPEDUM)
  EXPECT_EQ(2, dca::phys::domains::cluster_operations::origin_index(
                   unsorted_set_, dca::phys::domains::PARALLELLEPIPEDUM));

  // Sorted set (BRILLOUIN_ZONE)
  EXPECT_EQ(1, dca::phys::domains::cluster_operations::origin_index(
                   sorted_set_, dca::phys::domains::BRILLOUIN_ZONE));
}

TEST_F(ClusterOperationsTest, index) {
  // Unsorted set (PARALLELEPIPEDUM)
  for (int index = 0; index < unsorted_set_.size(); ++index) {
    EXPECT_EQ(index, dca::phys::domains::cluster_operations::index(
                         unsorted_set_[index], unsorted_set_, dca::phys::domains::PARALLELLEPIPEDUM));
  }

  // Sorted set (BRILLOUIN_ZONE)
  for (int index = 0; index < sorted_set_.size(); ++index) {
    EXPECT_EQ(index, dca::phys::domains::cluster_operations::index(
                         sorted_set_[index], sorted_set_, dca::phys::domains::BRILLOUIN_ZONE));
  }
}

#ifndef NDEBUG
TEST_F(ClusterOperationsDeathTest, index) {
  // Unsorted set (PARALLELEPIPEDUM)
  std::vector<double> not_element{-42., 24.};
  EXPECT_DEATH(dca::phys::domains::cluster_operations::index(not_element, unsorted_set_,
                                                             dca::phys::domains::PARALLELLEPIPEDUM),
               "index > -1 and index < elements.size()");

  // Sorted set (BRILLOUIN_ZONE)
  not_element = {-42., 24.};
  EXPECT_DEATH(dca::phys::domains::cluster_operations::index(not_element, sorted_set_,
                                                             dca::phys::domains::BRILLOUIN_ZONE),
               "math::util::distance2");
}
#endif  // NDEBUG

TEST_F(ClusterOperationsTest, translate_inside_cluster) {
  // Cluster vectors should be invariant.
  for (const auto& cluster_vector : r_cluster_elements_) {
    EXPECT_EQ(cluster_vector, dca::phys::domains::cluster_operations::translate_inside_cluster(
                                  cluster_vector, r_cluster_super_basis_));
  }

  // Translated cluster vector
  // {1., 0.} = {-1., 2.} + {2., -2.} = r_cluster_elements_[1] - r_super_basis_[1].
  std::vector<double> translated_vec{1., 0.};
  EXPECT_EQ(r_cluster_elements_[1], dca::phys::domains::cluster_operations::translate_inside_cluster(
                                        translated_vec, r_cluster_super_basis_));

  // Arbitrary vector inside the cluster should be invariant.
  std::vector<double> vec_inside{0., 0.5};
  EXPECT_EQ(vec_inside, dca::phys::domains::cluster_operations::translate_inside_cluster(
                            vec_inside, r_cluster_super_basis_));

  // Arbitrary vector outside the cluster.
  // {0., -0.5} = {0., 3.5} - {0., 4.} = {0., 3.5} - r_super_basis_[0] - r_super_basis_[1].
  std::vector<double> vec_outside{0., -0.5};
  std::vector<double> vec_outside_result{0., 3.5};
  EXPECT_EQ(vec_outside_result, dca::phys::domains::cluster_operations::translate_inside_cluster(
                                    vec_outside, r_cluster_super_basis_));
}

TEST_F(ClusterOperationsTest, find_closest_cluster_vector) {
  std::vector<double> q_input;
  std::vector<double> q_output;

  // Cluster vectors should be invariant.
  for (const auto& cluster_vector : r_cluster_elements_) {
    q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
        cluster_vector, r_cluster_elements_, r_cluster_super_basis_);
    EXPECT_EQ(cluster_vector, q_output);
  }

  // A translated cluster vector should return a cluster vector..
  q_input = {1., 0.};
  q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
      q_input, r_cluster_elements_, r_cluster_super_basis_);
  EXPECT_EQ(r_cluster_elements_[1], q_output);

  // Arbitray vector inside the cluster.
  q_input = {0., 0.25};
  q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
      q_input, r_cluster_elements_, r_cluster_super_basis_);
  EXPECT_EQ(r_cluster_elements_[2], q_output);

  // Arbitray vector outside the cluster.
  q_input = {0.75, 0.};
  q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
      q_input, r_cluster_elements_, r_cluster_super_basis_);
  EXPECT_EQ(r_cluster_elements_[1], q_output);

  // If there are multiple cluster vectors with the same minimal distance the first one in the
  // std::vector is returned.
  q_input = {0., 0.5};
  q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
      q_input, r_cluster_elements_, r_cluster_super_basis_);
  EXPECT_EQ(r_cluster_elements_[2], q_output);

  // Requiring a certain tolerance.
  double tolerance = 0.3;
  q_input = {0.75, 0.};
  EXPECT_NO_THROW(q_output = dca::phys::domains::cluster_operations::find_closest_cluster_vector(
                      q_input, r_cluster_elements_, r_cluster_super_basis_, tolerance));
  EXPECT_EQ(r_cluster_elements_[1], q_output);

  tolerance = 0.2;
  q_input = {0.75, 0.};
  EXPECT_THROW(dca::phys::domains::cluster_operations::find_closest_cluster_vector(
                   q_input, r_cluster_elements_, r_cluster_super_basis_, tolerance),
               std::logic_error);
}
