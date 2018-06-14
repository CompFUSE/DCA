// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements the cluster domain.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_HPP

#include <cassert>
#include <ios>
#include <vector>

#include "dca/linalg/matrix.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/cluster_specifications.hpp"
#include "dca/phys/domains/cluster/dual_cluster.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename cluster_type>
class cluster_symmetry {};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_domain {
public:
  const static int DIMENSION = D;
  using Scalar = scalar_type;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE SHAPE = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION =
      dual_cluster<REPRESENTATION>::REPRESENTATION;

  typedef cluster_domain<scalar_type, D, N, REPRESENTATION, S> this_type;
  typedef cluster_domain<scalar_type, D, N, DUAL_REPRESENTATION, S> dual_type;

  typedef
      typename cluster_specifications<scalar_type, N, R, S>::dmn_specifications_type dmn_specifications_type;

  typedef std::vector<scalar_type> element_type;

  static bool& is_initialized();

  static int& get_size();

  static std::vector<int>& get_dimensions();

  /*!
   *   The convention is that the basis and superbasis vectors are stored in the columns!
   */
  static scalar_type*& get_basis();
  static scalar_type*& get_super_basis();

  /*!
   *   The convention is that the inverse basis (and inverse superbasis) is defined as the inverse
   * of the basis (superbasis) matrix!
   */
  static scalar_type*& get_inverse_basis();
  static scalar_type*& get_inverse_super_basis();

  static std::vector<element_type>& get_basis_vectors();
  static std::vector<element_type>& get_super_basis_vectors();

  static std::string& get_name();

  static std::vector<element_type>& get_elements();

  static scalar_type& get_volume();

  static int origin_index();

  static int add(int i, int j);
  static int subtract(int i, int j);

  static dca::linalg::Matrix<int, dca::linalg::CPU>& get_add_matrix();
  static dca::linalg::Matrix<int, dca::linalg::CPU>& get_subtract_matrix();

  static void reset();

  template <typename ss_type>
  static void print(ss_type& ss);
};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool& cluster_domain<scalar_type, D, N, R, S>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int& cluster_domain<scalar_type, D, N, R, S>::get_size() {
  static int size = 0;
  return size;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<int>& cluster_domain<scalar_type, D, N, R, S>::get_dimensions() {
  static std::vector<int> dimensions;
  return dimensions;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, D, N, R, S>::get_basis() {
  static scalar_type* basis = NULL;
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, D, N, R, S>::get_super_basis() {
  static scalar_type* basis = NULL;
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, D, N, R, S>::get_inverse_basis() {
  static scalar_type* basis = NULL;
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& cluster_domain<scalar_type, D, N, R, S>::get_inverse_super_basis() {
  static scalar_type* basis = NULL;
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& cluster_domain<scalar_type, D, N, R, S>::get_basis_vectors() {
  static std::vector<element_type> basis;
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& cluster_domain<scalar_type, D, N, R, S>::get_super_basis_vectors() {
  static std::vector<element_type> super_basis;
  return super_basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::string& cluster_domain<scalar_type, D, N, R, S>::get_name() {
  static std::string name = to_str(NAME) + " " + to_str(REPRESENTATION) + " " + to_str(SHAPE) +
                            " (DIMENSION : " + to_str(DIMENSION) + ")";
  return name;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& cluster_domain<scalar_type, D, N, R, S>::get_elements() {
  static std::vector<element_type> elements(get_size());
  return elements;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type& cluster_domain<scalar_type, D, N, R, S>::get_volume() {
  static scalar_type volume = 0;
  return volume;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, D, N, R, S>::origin_index() {
  static int index = cluster_operations::origin_index(get_elements(), SHAPE);
  assert(math::util::l2Norm2(get_elements()[index]) < 1.e-6);
  return index;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, D, N, R, S>::add(int i, int j) {
  static dca::linalg::Matrix<int, dca::linalg::CPU>& A = get_add_matrix();
  return A(i, j);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
dca::linalg::Matrix<int, dca::linalg::CPU>& cluster_domain<scalar_type, D, N, R, S>::get_add_matrix() {
  assert(SHAPE == BRILLOUIN_ZONE);
  static dca::linalg::Matrix<int, dca::linalg::CPU> A("add", get_size());
  return A;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int cluster_domain<scalar_type, D, N, R, S>::subtract(int i, int j) {
  static dca::linalg::Matrix<int, dca::linalg::CPU>& A = get_subtract_matrix();
  return A(i, j);
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
dca::linalg::Matrix<int, dca::linalg::CPU>& cluster_domain<scalar_type, D, N, R, S>::get_subtract_matrix() {
  assert(SHAPE == BRILLOUIN_ZONE);
  static dca::linalg::Matrix<int, dca::linalg::CPU> A("subtract", get_size());
  return A;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void cluster_domain<scalar_type, D, N, R, S>::reset() {
  get_size() = 0;
  get_name() = "";

  get_elements().resize(0);

  if (get_basis() != NULL)
    delete[] get_basis();

  if (get_super_basis() != NULL)
    delete[] get_super_basis();

  if (get_inverse_basis() != NULL)
    delete[] get_inverse_basis();

  if (get_inverse_super_basis() != NULL)
    delete[] get_inverse_super_basis();

  is_initialized() = false;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
template <typename ss_type>
void cluster_domain<scalar_type, D, N, R, S>::print(ss_type& ss) {
  ss << std::scientific;
  ss.precision(6);

  ss << "\t name        : " << get_name() << "\n";
  ss << "\t name (dual) : " << dual_type::get_name() << "\n\n";

  ss << "\t size        : " << get_size() << "\n\n";

  ss << "\t\t\t" << to_str(REPRESENTATION) << "\t\t\t|\t" << to_str(DUAL_REPRESENTATION) << "\n";
  ss << "\t origin-index : " << origin_index() << "\t\t\t\t|\t" << dual_type::origin_index() << "\n";
  ss << "\t volume       : " << get_volume() << "\t\t\t|\t" << dual_type::get_volume() << "\n\n";

  ss << "\t basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "|\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << dual_type::get_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t super-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "|\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << dual_type::get_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t inverse-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_inverse_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "|\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << dual_type::get_inverse_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t inverse-super-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_inverse_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "|\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << dual_type::get_inverse_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n\n";

  if (NAME == CLUSTER) {
    for (int l = 0; l < get_size(); l++) {
      ss << "\t" << l << "\t|\t";
      math::util::print(get_elements()[l]);
      ss << "\t";
      math::util::print(dual_type::get_elements()[l]);
      ss << "\n";
    }
    ss << "\n\n\t" << to_str(REPRESENTATION) << " k-space symmetries : \n\n";
    {
      for (int i = 0; i < this_type::get_size(); ++i) {
        for (int j = 0; j < cluster_symmetry<this_type>::b_dmn_t::dmn_size(); ++j) {
          ss << "\t" << i << ", " << j << "\t|\t";

          for (int l = 0; l < cluster_symmetry<this_type>::sym_super_cell_dmn_t::dmn_size(); ++l)
            ss << "\t" << cluster_symmetry<this_type>::get_symmetry_matrix()(i, j, l).first << ", "
               << cluster_symmetry<this_type>::get_symmetry_matrix()(i, j, l).second;

          ss << "\n";
        }
      }
      ss << "\n";
    }

    ss << "\n\n\t" << to_str(DUAL_REPRESENTATION) << " symmetries : \n\n";
    {
      for (int i = 0; i < dual_type::get_size(); ++i) {
        for (int j = 0; j < cluster_symmetry<dual_type>::b_dmn_t::dmn_size(); ++j) {
          ss << "\t" << i << ", " << j << "\t|\t";

          for (int l = 0; l < cluster_symmetry<dual_type>::sym_super_cell_dmn_t::dmn_size(); ++l)
            ss << "\t" << cluster_symmetry<dual_type>::get_symmetry_matrix()(i, j, l).first << ", "
               << cluster_symmetry<dual_type>::get_symmetry_matrix()(i, j, l).second;

          ss << "\n";
        }
      }
      ss << "\n";
    }
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_HPP
