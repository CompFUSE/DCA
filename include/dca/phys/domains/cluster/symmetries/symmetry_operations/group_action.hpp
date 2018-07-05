// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Group action.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_GROUP_ACTION_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_GROUP_ACTION_HPP

#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int DIMENSION>
class group_action {
  typedef group_action<DIMENSION> this_type;

public:
  // The transformation-matrix is stored here as a linear array in column major order.

  // action : transforms vector vec -> symmetry::matrix().vec
  static void action(double* r_vec, const double* matrix);
  static void action(std::vector<double>& vec, const double* matrix);

  //  vec1 == symmetry::matrix().vec2 ?
  static bool is_equivalent(double* vec1, double* vec2, const double* matrix);
  static bool is_equivalent(std::vector<double>& vec1, std::vector<double>& vec2,
                            const double* matrix);
};

template <int DIMENSION>
void group_action<DIMENSION>::action(double* /*vec*/, const double* /*matrix*/) {}

template <>
void group_action<2>::action(double* vec, const double* matrix);
template <>
void group_action<3>::action(double* vec, const double* matrix);

template <>
void group_action<2>::action(std::vector<double>& vec, const double* matrix);
template <>
void group_action<3>::action(std::vector<double>& vec, const double* matrix);

template <int DIMENSION>
bool group_action<DIMENSION>::is_equivalent(double* /*vec1*/, double* /*vec2*/,
                                            const double* /*matrix*/) {
  return false;
}

template <>
bool group_action<2>::is_equivalent(double* vec1, double* vec2, const double* matrix);
template <>
bool group_action<3>::is_equivalent(double* vec1, double* vec2, const double* matrix);

template <>
bool group_action<2>::is_equivalent(std::vector<double>& vec1, std::vector<double>& vec2,
                                    const double* matrix);
template <>
bool group_action<3>::is_equivalent(std::vector<double>& vec1, std::vector<double>& vec2,
                                    const double* matrix);

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_GROUP_ACTION_HPP
