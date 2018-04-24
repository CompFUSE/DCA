// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This lattice allows to simulate a 1D chain by setting the hopping in one dimension to zero.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SINGLEBAND_CHAIN_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SINGLEBAND_CHAIN_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class singleband_chain {
public:
  typedef domains::no_symmetry<1> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 1;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are two-dimensional (access through index 0 and 1).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
};

template <typename point_group_type>
double* singleband_chain<point_group_type>::initialize_r_DCA_basis() {
  static std::array<double, 4> r_DCA{1., 0., 0., 1};
  return r_DCA.data();
}

template <typename point_group_type>
double* singleband_chain<point_group_type>::initialize_k_DCA_basis() {
  static std::array<double, 4> k_DCA{2. * M_PI, 0., 0, 2 * M_PI};
  return k_DCA.data();
}

template <typename point_group_type>
double* singleband_chain<point_group_type>::initialize_r_LDA_basis() {
  static std::array<double, 4> r_LDA{1., 0., 0., 1};

  return r_LDA.data();
}

template <typename point_group_type>
double* singleband_chain<point_group_type>::initialize_k_LDA_basis() {
  static std::array<double, 4> k_LDA{2. * M_PI, 0., 0, 2 * M_PI};
  return k_LDA.data();
}

template <typename point_group_type>
std::vector<int> singleband_chain<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> singleband_chain<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs{std::vector<double>{0,0}};

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> singleband_chain<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void singleband_chain<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Singleband chainhas one band.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // Same band, opposite spin.
  H_interaction = 0.;

  for (int b1 = 0; b1 < BANDS; b1++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int b2 = 0; b2 < BANDS; b2++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (b1 == b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;
        }
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void singleband_chain<point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void singleband_chain<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("lattice has one bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t = parameters.get_t();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto k = k_vecs[k_ind];
    const auto val = -2. * t[0] * std::cos(k[0]) - 2. * t[1] * std::cos(k[1]);

    for (int s = 0; s < 2; ++s)
      H_0(0, s, 0, s, k_ind) = val;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SINGLEBAND_CHAIN_HPP
