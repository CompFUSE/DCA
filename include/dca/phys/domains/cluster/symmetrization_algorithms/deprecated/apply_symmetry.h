// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Apply symmetry.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRY_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRY_H

#include <cmath>
#include <vector>

#include "dca/util/type_list.hpp"

template <class cluster_type, class point_group, int INDEX>
struct apply_symmetry {
  const static int DIMENSION = cluster_type::DIMENSION;

  typedef typename cluster_type::base_cluster base_cluster_type;

  typedef typename dca::util::TypeAt<INDEX, point_group>::Result symmetry_type;
  typedef typename symmetry_type::base_type group_action_type;

  template <class cluster_reduction_type>
  static void reduce_cluster_coordinates(cluster_reduction_type& cluster_reduction) {
    bool is_invariant = cluster_is_invariant_under_symmetry();

    if (is_invariant) {
      for (int i = cluster_type::get_size() - 1; i >= 0; i--) {
        if (i == cluster_reduction.get_irreducible_index(i)) {
          std::vector<double> vec_1 = cluster_type::get_elements()[i];

          group_action_type::action(vec_1, symmetry_type::matrix());

          std::vector<double> vec_2 = cluster_type::back_inside_cluster(vec_1);

          int index = search(vec_2, cluster_reduction);

          while (index != cluster_reduction.get_irreducible_index(index))
            index = cluster_reduction.get_irreducible_index(index);

          if (index != i) {
            cluster_reduction.get_irreducible_index(i) = index;

            cluster_reduction.get_weight(index) += cluster_reduction.get_weight(i);
            cluster_reduction.get_weight(i) = 0;
          }
        }
      }

      reduce(cluster_reduction);
    }

    apply_symmetry<cluster_type, point_group, (INDEX - 1)>::reduce_cluster_coordinates(
        cluster_reduction);
  }

  template <class cluster_reduction_type>
  static void reduce_cluster_pair(int i, cluster_reduction_type& cluster_reduction) {
    bool is_invariant = cluster_is_invariant_under_symmetry();

    if (is_invariant) {
      std::vector<double> vec_1 = cluster_type::get_elements()[i];

      group_action_type::action(vec_1, symmetry_type::matrix());

      std::vector<double> vec_2 = cluster_type::back_inside_cluster(vec_1);

      int index = search(vec_2, cluster_reduction);

      if (index == cluster_reduction.get_irreducible_index(i))
        apply_symmetry<cluster_type, point_group, INDEX>::reduce_pair_with_this_symmetry(
            i, cluster_reduction);
    }

    apply_symmetry<cluster_type, point_group, INDEX - 1>::reduce_cluster_pair(i, cluster_reduction);
  }

  template <class cluster_reduction_type>
  static void reduce_pair_with_this_symmetry(int i, cluster_reduction_type& cluster_reduction) {
    std::vector<double> vec_i_1 = cluster_type::get_elements()[i];
    group_action_type::action(vec_i_1, symmetry_type::matrix());
    std::vector<double> vec_i_2 = cluster_type::back_inside_cluster(vec_i_1);
    int index_i = search(vec_i_2, cluster_reduction);

    for (int j = 0; j < cluster_type::get_size(); j++) {
      cluster_reduction.get_irreducible_pair(i, j).first = index_i;

      std::vector<double> vec_j_1 = cluster_type::get_elements()[j];
      group_action_type::action(vec_j_1, symmetry_type::matrix());
      std::vector<double> vec_j_2 = cluster_type::back_inside_cluster(vec_j_1);
      int index_j = search(vec_j_2, cluster_reduction);

      cluster_reduction.get_irreducible_pair(i, j).second = index_j;
    }
  }

  static bool cluster_is_invariant_under_symmetry() {
    // if the cluster is invariant under the symmetry operation O,
    // then the action on the spanning vectors should only be a
    // permutation of the spanning vectors...
    // --> Action[ super_basis_vectors] = Permutation * super_basis_vectors

    double* permutation = new double[DIMENSION * DIMENSION];

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        permutation[i + j * DIMENSION] = 0;

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        for (int l1 = 0; l1 < DIMENSION; l1++)
          for (int l2 = 0; l2 < DIMENSION; l2++)
            permutation[i + j * DIMENSION] += base_cluster_type::get_k_super_basis()[i][l1] *
                                              symmetry_type::matrix()[l1 + DIMENSION * l2] *
                                              base_cluster_type::get_r_super_basis()[j][l2] /
                                              (2. * M_PI);

    bool is_permutation_matrix = true;

    for (int i = 0; i < DIMENSION; i++)
      for (int j = 0; j < DIMENSION; j++)
        if (!(std::fabs(std::fabs(permutation[i + j * DIMENSION]) - 1) < 1.e-6 ||
              std::fabs(permutation[i + j * DIMENSION]) < 1.e-6))
          is_permutation_matrix = false;

    for (int i = 0; i < DIMENSION; i++) {
      double result_col = 0;
      double result_row = 0;

      for (int j = 0; j < DIMENSION; j++) {
        result_col += std::fabs(permutation[j + i * DIMENSION]);
        result_row += std::fabs(permutation[i + j * DIMENSION]);
      }

      if (!(std::fabs(result_col - 1) < 1.e-6 && std::fabs(result_col - 1) < 1.e-6))
        is_permutation_matrix = false;
    }

    delete[] permutation;

    return is_permutation_matrix;
  }

  template <class cluster_reduction_type>
  static int search(std::vector<double>& vec, cluster_reduction_type& cluster_reduction) {
    for (int i = 0; i < cluster_type::get_size(); i++) {
      if (L2_norm(vec, cluster_type::get_elements()[i]) < 1.e-6)
        return i;
    }

    throw std::logic_error(__FUNCTION__);
  }

  template <class cluster_reduction_type>
  static void reduce(cluster_reduction_type& cluster_reduction) {
    for (int i = 0; i < cluster_type::get_size(); i++) {
      int index = i;

      while (index != cluster_reduction.get_irreducible_index(index))
        index = cluster_reduction.get_irreducible_index(index);

      cluster_reduction.get_irreducible_index(i) = index;
    }
  }
};

/*!
 *  \ingroup SYMMETRIES
 *
 *  \author  Peter Staar
 */
template <class cluster_type, class point_group>
struct apply_symmetry<cluster_type, point_group, -1> {
  template <class cluster_reduction_type>
  static void reduce_cluster_coordinates(cluster_reduction_type& cluster_reduction) {}

  template <class cluster_reduction_type>
  static void reduce_cluster_pair(int i, cluster_reduction_type& cluster_reduction) {}
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRY_H
