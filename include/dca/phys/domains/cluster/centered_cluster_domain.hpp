// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements the centered cluster domain.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CENTERED_CLUSTER_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CENTERED_CLUSTER_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename cluster_type>
class centered_cluster_domain {};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>> {
public:
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE SHAPE = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION =
      dual_cluster<REPRESENTATION>::REPRESENTATION;

  typedef cluster_domain<scalar_type, D, N, REPRESENTATION, S> cluster_type;

  typedef centered_cluster_domain<cluster_domain<scalar_type, D, N, REPRESENTATION, S>> this_type;
  typedef centered_cluster_domain<cluster_domain<scalar_type, D, N, DUAL_REPRESENTATION, S>> dual_type;

  typedef
      typename cluster_specifications<scalar_type, N, R, S>::dmn_specifications_type dmn_specifications_type;

  typedef std::vector<scalar_type> element_type;

  static bool& is_initialized();

  static int& get_size();

  static std::vector<int>& get_dimensions();

  static scalar_type*& get_basis();
  static scalar_type*& get_super_basis();

  static scalar_type*& get_inverse_basis();
  static scalar_type*& get_inverse_super_basis();

  static std::vector<element_type>& get_basis_vectors();
  static std::vector<element_type>& get_super_basis_vectors();

  static std::string& get_name();

  static std::vector<scalar_type>& get_weights();
  static std::vector<element_type>& get_elements();

  static scalar_type& get_volume();

  static int origin_index();

  static void initialize();

  static void reset();

  template <typename ss_type>
  static void print(ss_type& ss);
};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
bool& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_size() {
  static int size = 0;
  return size;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<int>& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_dimensions() {
  static std::vector<int> dimensions = cluster_type::get_dimensions();
  return dimensions;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_basis() {
  static scalar_type* basis = cluster_type::get_basis();
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_super_basis() {
  static scalar_type* basis = cluster_type::get_super_basis();
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_inverse_basis() {
  static scalar_type* basis = cluster_type::get_inverse_basis();
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type*& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_inverse_super_basis() {
  static scalar_type* basis = cluster_type::get_inverse_super_basis();
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& centered_cluster_domain<
    cluster_domain<scalar_type, D, N, R, S>>::get_basis_vectors() {
  static std::vector<element_type> basis = cluster_type::get_basis_vectors();
  return basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& centered_cluster_domain<
    cluster_domain<scalar_type, D, N, R, S>>::get_super_basis_vectors() {
  static std::vector<element_type> super_basis = cluster_type::get_super_basis_vectors();
  return super_basis;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::string& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_name() {
  static std::string name = "centered-cluster-domain " + to_str(NAME) + " " + to_str(REPRESENTATION) +
                            " " + to_str(SHAPE) + " (DIMENSION : " + to_str(DIMENSION) + ")";
  return name;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<scalar_type>& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_weights() {
  static std::vector<scalar_type> weights(get_size());
  return weights;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
std::vector<std::vector<scalar_type>>& centered_cluster_domain<
    cluster_domain<scalar_type, D, N, R, S>>::get_elements() {
  static std::vector<element_type> elements(get_size());
  return elements;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
scalar_type& centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::get_volume() {
  static scalar_type volume = cluster_type::get_volume();
  return volume;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
int centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::origin_index() {
  static int index = cluster_operations::origin_index(get_elements());
  assert(math::util::l2Norm2(get_elements()[index]) < 1.e-6);
  return index;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::initialize() {
  if (not is_initialized()) {
    get_weights().resize(0);
    get_elements().resize(0);

    for (int i = 0; i < cluster_type::get_size(); i++) {
      std::vector<double> vec = cluster_type::get_elements()[i];

      std::vector<std::vector<double>> vecs =
          cluster_operations::equivalent_vectors(vec, cluster_type::get_super_basis_vectors());

      for (size_t j = 0; j < vecs.size(); j++) {
        get_weights().push_back(1. / vecs.size());
        get_elements().push_back(vecs[j]);
      }

      // 	get_weights() .push_back(1.);
      // 	get_elements().push_back(vecs[0]);
    }

    get_size() = get_elements().size();

    is_initialized() = true;
  }
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
void centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::reset() {
  get_size() = 0;
  get_name() = "";

  get_elements().resize(0);

  is_initialized() = false;
}

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
template <typename ss_type>
void centered_cluster_domain<cluster_domain<scalar_type, D, N, R, S>>::print(ss_type& ss) {
  ss << std::scientific;
  ss.precision(6);

  ss << "\t name        : " << get_name() << "\n";

  ss << "\t size        : " << get_size() << "\n\n";

  // ss << "\t origin-index : " << origin_index() << "\n";
  ss << "\t volume       : " << get_volume() << "\n\n";

  ss << "\t basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t super-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t inverse-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_inverse_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n";

  ss << "\t inverse-super-basis : \n";
  for (int d0 = 0; d0 < DIMENSION; d0++) {
    ss << "\t\t\t";
    for (int d1 = 0; d1 < DIMENSION; d1++)
      ss << get_inverse_super_basis()[d0 + d1 * DIMENSION] << "\t";
    ss << "\n";
  }
  ss << "\n\n";

  for (int l = 0; l < get_size(); l++) {
    ss << "\t" << l << "\t";
    math::util::print(get_elements()[l]);
    ss << "\t" << get_weights()[l] << "\n";
  }
  ss << "\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CENTERED_CLUSTER_DOMAIN_HPP
