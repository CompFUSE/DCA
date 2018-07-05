// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a cluster domain iterator.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_ITERATOR_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_ITERATOR_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename cluster_type>
class cluster_domain_iterator {};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>> {
public:
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE SHAPE = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION =
      dual_cluster<REPRESENTATION>::REPRESENTATION;

  typedef cluster_domain<scalar_type, D, N, REAL_SPACE, S> r_dmn;
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_dmn;

  typedef cluster_domain<scalar_type, D, TMP_CLUSTER, REAL_SPACE, S> tmp_r_dmn;
  typedef cluster_domain<scalar_type, D, TMP_CLUSTER, MOMENTUM_SPACE, S> tmp_k_dmn;

  typedef std::vector<scalar_type> element_type;

  cluster_domain_iterator();

  void execute(int N_MAX, double l_min, double theta_min);

private:
  void execute_2D(int N_MAX, double l_min, double theta_min);

  bool find_anti_nodal_points();

  bool find_pi_zero();
  bool find_zero_pi();

  struct cluster_info_struct {
    int size;

    double theta;

    double l0;
    double l1;

    bool pi_zero;
    bool zero_pi;

    std::vector<int> vec_0;
    std::vector<int> vec_1;

    static bool comparison(cluster_info_struct i, cluster_info_struct j) {
      if (i.size == j.size) {
        if (i.pi_zero == j.pi_zero) {
          if (i.zero_pi == j.zero_pi) {
            if (i.vec_0 == j.vec_0)
              return math::util::isLessVector(i.vec_1, j.vec_1);
            else
              return math::util::isLessVector(i.vec_0, j.vec_0);
          }
          else
            return i.zero_pi > j.zero_pi;
        }
        else
          return i.pi_zero > j.pi_zero;
      }
      else
        return (i.size < j.size);
    }
  };

  std::vector<cluster_info_struct> list_of_clusters;
};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::cluster_domain_iterator()
    : list_of_clusters(0) {}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::execute(int N_MAX,
                                                                               double l_min,
                                                                               double theta_min) {
  switch (DIMENSION) {
    //     case 1:
    //       execute_1D();
    //       break;

    case 2:
      execute_2D(N_MAX, l_min, theta_min);
      break;

    //     case 3:
    //       execute_3D(shift);
    //       break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::execute_2D(int N_MAX,
                                                                                  double l_min,
                                                                                  double theta_min) {
  for (int d0_x = 0; d0_x <= N_MAX; d0_x++) {
    for (int d0_y = 0; d0_y <= d0_x; d0_y++) {
      for (int d1_x = -N_MAX; d1_x <= N_MAX; d1_x++) {
        for (int d1_y = 0; d1_y <= N_MAX; d1_y++) {
          std::vector<int> vec_0(DIMENSION, 0);
          std::vector<int> vec_1(DIMENSION, 0);

          vec_0[0] = d0_x;
          vec_0[1] = d0_y;

          vec_1[0] = d1_x;
          vec_1[1] = d1_y;

          if (abs(vec_0[0] * vec_1[1] - vec_0[1] * vec_1[0]) > 1.e-6) {
            tmp_r_dmn::reset();

            std::vector<std::vector<int>> basis;
            basis.push_back(vec_0);
            basis.push_back(vec_1);

            cluster_domain_initializer<func::dmn_0<tmp_r_dmn>>::execute(r_dmn::get_basis(), basis,
                                                                        false);

            std::vector<double> b_0 = tmp_r_dmn::get_super_basis_vectors()[0];
            std::vector<double> b_1 = tmp_r_dmn::get_super_basis_vectors()[1];

            double l0 = std::sqrt(b_0[0] * b_0[0] + b_0[1] * b_0[1]) + 1.e-6;
            double l1 = std::sqrt(b_1[0] * b_1[0] + b_1[1] * b_1[1]) + 1.e-6;

            double theta =
                std::acos(std::abs(b_0[0] * b_1[0] + b_0[1] * b_1[1]) / (l0 * l1)) * 180.0 / M_PI;

            if (find_anti_nodal_points() and std::abs(l0 / l1) > 1. - l_min and
                std::abs(l0 / l1) < 1. + l_min and theta > 90 - theta_min and theta < 90 + theta_min and
                find_pi_zero() and (tmp_r_dmn::get_size() == 20 or tmp_r_dmn::get_size() == 24 or
                                    tmp_r_dmn::get_size() == 28 or tmp_r_dmn::get_size() == 32 or
                                    tmp_r_dmn::get_size() == 36 or tmp_r_dmn::get_size() == 40 or
                                    tmp_r_dmn::get_size() == 44 or tmp_r_dmn::get_size() == 48 or
                                    tmp_r_dmn::get_size() == 52 or tmp_r_dmn::get_size() == 56)) {
              cluster_info_struct cluster;
              cluster.size = tmp_r_dmn::get_size();

              cluster.l0 = l0;
              cluster.l1 = l1;

              cluster.theta = theta;

              cluster.pi_zero = find_pi_zero();
              cluster.zero_pi = find_zero_pi();

              cluster.vec_0 = vec_0;
              cluster.vec_1 = vec_1;

              list_of_clusters.push_back(cluster);
            }
          }
        }
      }
    }
  }

  std::sort(list_of_clusters.begin(), list_of_clusters.end(), cluster_info_struct::comparison);

  for (size_t l = 0; l < list_of_clusters.size(); l++) {
    std::cout << "[" << list_of_clusters[l].size << ",\t";  //<< list_of_clusters[l].l0 << "\t" <<
                                                            // list_of_clusters[l].l1 << "\t" <<
    // list_of_clusters[l].theta << "\t";
    std::cout << list_of_clusters[l].pi_zero << ",\t" << list_of_clusters[l].zero_pi << ",\t";
    std::cout << "[" << list_of_clusters[l].vec_0[0] << ", " << list_of_clusters[l].vec_0[1]
              << "],  ";
    std::cout << "[" << list_of_clusters[l].vec_1[0] << ", " << list_of_clusters[l].vec_1[1]
              << "] ";
    std::cout << "],\n";
  }
  std::cout << "\n";
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::find_anti_nodal_points() {
  bool contains = false;

  std::vector<double> pi_zero(2, 0);
  std::vector<double> zero_pi(2, 0);

  pi_zero[0] = M_PI;
  zero_pi[1] = M_PI;

  for (int l = 0; l < tmp_k_dmn::get_size(); l++) {
    double dif0 = math::util::distance(pi_zero, tmp_k_dmn::get_elements()[l]);
    double dif1 = math::util::distance(zero_pi, tmp_k_dmn::get_elements()[l]);

    if (dif0 < 1.e-6 or dif1 < 1.e-6)
      contains = true;
  }

  return contains;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::find_pi_zero() {
  bool contains = false;

  std::vector<double> pi_zero(2, 0);

  pi_zero[0] = M_PI;

  for (int l = 0; l < tmp_k_dmn::get_size(); l++) {
    double dif0 = math::util::distance(pi_zero, tmp_k_dmn::get_elements()[l]);

    if (dif0 < 1.e-6)
      contains = true;
  }

  return contains;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool cluster_domain_iterator<cluster_domain<scalar_type, D, N, R, S>>::find_zero_pi() {
  bool contains = false;

  std::vector<double> zero_pi(2, 0);

  zero_pi[1] = M_PI;

  for (int l = 0; l < tmp_k_dmn::get_size(); l++) {
    double dif1 = math::util::distance(zero_pi, tmp_k_dmn::get_elements()[l]);

    if (dif1 < 1.e-6)
      contains = true;
  }

  return contains;
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_ITERATOR_HPP
