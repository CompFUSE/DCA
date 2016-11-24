// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Bilayer lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_LATTICE_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class bilayer_lattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 2;

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
  template <typename ParametersType, typename ScalarType, typename BandSpinDmn, typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KDmn>>& H_0);
};

template <typename point_group_type>
double* bilayer_lattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* bilayer_lattice<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* bilayer_lattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* bilayer_lattice<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> bilayer_lattice<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> bilayer_lattice<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> bilayer_lattice<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void bilayer_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();              // Same band, opposite spin.
  const double V = parameters.get_V();              // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();  // Different band, same spin.

  H_interaction = 0.;

  for (int b1 = 0; b1 < BANDS; b1++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int b2 = 0; b2 < BANDS; b2++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (b1 == b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;

          if (b1 != b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = V;

          if (b1 != b2 && s1 == s2)
            H_interaction(b1, s1, b2, s2, origin) = V_prime;
        }
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void bilayer_lattice<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandSpinDmn, typename KDmn>
void bilayer_lattice<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KDmn>>& H_0) {
  static_assert(util::Length<typename BandSpinDmn::this_type>::value == 2,
                "BandSpinDmn has two subdomains, the band domain and the spin domain.");

  using BandDmn = typename util::TypeAt<0, typename BandSpinDmn::template domain_typelist<0>>::type;
  using SpinDmn = typename util::TypeAt<0, typename BandSpinDmn::template domain_typelist<1>>::type;

  if (BandDmn::get_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::get_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t = parameters.get_t();
  const auto t_prime = parameters.get_t_prime();
  const auto t_perp = parameters.get_t_perp();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const auto val =
        -2. * t * (std::cos(k[0]) + std::cos(k[1])) - 4. * t_prime * std::cos(k[0]) * std::cos(k[1]);

    H_0(0, 0, 0, 0, k_ind) = val;
    H_0(0, 1, 0, 1, k_ind) = val;
    H_0(1, 0, 1, 0, k_ind) = val;
    H_0(1, 1, 1, 1, k_ind) = val;

    H_0(0, 0, 1, 0, k_ind) = -t_perp;
    H_0(0, 1, 1, 1, k_ind) = -t_perp;
    H_0(1, 0, 0, 0, k_ind) = -t_perp;
    H_0(1, 1, 0, 1, k_ind) = -t_perp;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_LATTICE_HPP
