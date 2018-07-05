// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements group_action.hpp.

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"
#include <cmath>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <>
void group_action<2>::action(double* vec, const double* matrix) {
  double tmp_vec[2] = {0., 0.};

  tmp_vec[0] = matrix[0] * vec[0] + matrix[2] * vec[1];
  tmp_vec[1] = matrix[1] * vec[0] + matrix[3] * vec[1];

  vec[0] = tmp_vec[0];
  vec[1] = tmp_vec[1];
}

template <>
void group_action<3>::action(double* vec, const double* matrix) {
  double tmp_vec[3] = {0., 0., 0.};

  tmp_vec[0] = matrix[0] * vec[0] + matrix[3] * vec[1] + matrix[6] * vec[2];
  tmp_vec[1] = matrix[1] * vec[0] + matrix[4] * vec[1] + matrix[7] * vec[2];
  tmp_vec[2] = matrix[2] * vec[0] + matrix[5] * vec[1] + matrix[8] * vec[2];

  vec[0] = tmp_vec[0];
  vec[1] = tmp_vec[1];
  vec[2] = tmp_vec[2];
}

template <>
void group_action<2>::action(std::vector<double>& vec, const double* matrix) {
  double tmp_vec[2] = {0., 0.};

  tmp_vec[0] = matrix[0] * vec[0] + matrix[2] * vec[1];
  tmp_vec[1] = matrix[1] * vec[0] + matrix[3] * vec[1];

  vec[0] = tmp_vec[0];
  vec[1] = tmp_vec[1];
}

template <>
void group_action<3>::action(std::vector<double>& vec, const double* matrix) {
  double tmp_vec[3] = {0., 0., 0.};

  tmp_vec[0] = matrix[0] * vec[0] + matrix[3] * vec[1] + matrix[6] * vec[2];
  tmp_vec[1] = matrix[1] * vec[0] + matrix[4] * vec[1] + matrix[7] * vec[2];
  tmp_vec[2] = matrix[2] * vec[0] + matrix[5] * vec[1] + matrix[8] * vec[2];

  vec[0] = tmp_vec[0];
  vec[1] = tmp_vec[1];
  vec[2] = tmp_vec[2];
}

template <>
bool group_action<2>::is_equivalent(double* vec1, double* vec2, const double* matrix) {
  group_action<2>::action(vec2, matrix);

  if ((std::pow(vec1[0] - vec2[0], 2.) + std::pow(vec1[1] - vec2[1], 2.)) < 1.e-6)
    return true;
  else
    return false;
}

template <>
bool group_action<3>::is_equivalent(double* vec1, double* vec2, const double* matrix) {
  group_action<3>::action(vec2, matrix);

  if ((std::pow(vec1[0] - vec2[0], 2.) + std::pow(vec1[1] - vec2[1], 2.) +
       std::pow(vec1[2] - vec2[2], 2.)) < 1.e-6)
    return true;
  else
    return false;
}

template <>
bool group_action<2>::is_equivalent(std::vector<double>& vec1, std::vector<double>& vec2,
                                    const double* matrix) {
  std::vector<double> transformed_vec = vec2;
  group_action<2>::action(transformed_vec, matrix);

  if ((std::pow(vec1[0] - transformed_vec[0], 2.) + std::pow(vec1[1] - transformed_vec[1], 2.)) <
      1.e-6)
    return true;
  else
    return false;
}

template <>
bool group_action<3>::is_equivalent(std::vector<double>& vec1, std::vector<double>& vec2,
                                    const double* matrix) {
  std::vector<double> transformed_vec = vec2;
  group_action<3>::action(transformed_vec, matrix);

  if ((std::pow(vec1[0] - transformed_vec[0], 2.) + std::pow(vec1[1] - transformed_vec[1], 2.) +
       std::pow(vec1[2] - transformed_vec[2], 2.)) < 1.e-6)
    return true;
  else
    return false;
}

}  // domains
}  // phys
}  // dca
