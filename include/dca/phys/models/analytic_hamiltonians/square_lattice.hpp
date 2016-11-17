// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Square lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP

#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class square_lattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static int DIMENSION = 2;
  const static int BANDS = 1;

  static double* initialize_r_DCA_basis();
  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

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
double* square_lattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.;
  r_DCA[1] = 0.;
  r_DCA[2] = 0.;
  r_DCA[3] = 1.;

  return r_DCA;
}

template <typename point_group_type>
double* square_lattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
std::vector<int> square_lattice<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> square_lattice<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> square_lattice<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

// TODO: Add non-local interaction of same spins.
//       Use V instead of U_prime?
template <typename point_group_type>
template <class domain, class parameters_type>
void square_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_interaction, parameters_type& parameters) {
  H_interaction = 0.;

  // actually the same as DCA_r_cluster_type (see typedifinitions.h).
  typedef
      typename dca::util::TypeAt<0, typename domain::template domain_typelist<2>>::type DCA_r_cluster_t;

  int DIMENSION = DCA_r_cluster_t::DIMENSION;
  assert(DIMENSION == 2);

  int origin = DCA_r_cluster_t::origin_index();

  std::vector<typename DCA_r_cluster_t::element_type>& basis = DCA_r_cluster_t::get_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& super_basis =
      DCA_r_cluster_t::get_super_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& elements = DCA_r_cluster_t::get_elements();

  std::vector<int> nn_index(DIMENSION);  // Indices of nearest neighbours w.r.t. origin.
  for (int d = 0; d < DIMENSION; ++d) {
    std::vector<double> basis_vec =
        domains::cluster_operations::translate_inside_cluster(basis[d], super_basis);
    nn_index[d] = domains::cluster_operations::index(basis_vec, elements, domains::BRILLOUIN_ZONE);
  }

  // Nearest-neighbor opposite spin interaction
  double V = parameters.get_V();
  H_interaction(0, 1, nn_index[0]) = V;
  H_interaction(1, 0, nn_index[0]) = V;
  H_interaction(0, 1, nn_index[1]) = V;
  H_interaction(1, 0, nn_index[1]) = V;

  // Nearest-neighbor same spin interaction
  double V_prime = parameters.get_V_prime();
  H_interaction(0, 0, nn_index[0]) = V_prime;
  H_interaction(1, 1, nn_index[0]) = V_prime;
  H_interaction(0, 0, nn_index[1]) = V_prime;
  H_interaction(1, 1, nn_index[1]) = V_prime;

  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  double U = parameters.get_U();
  H_interaction(0, 1, origin) = U;
  H_interaction(1, 0, origin) = U;
}

template <typename point_group_type>
template <class domain>
void square_lattice<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries(0, 0) = 0;
  H_symmetries(0, 1) = -1;
  H_symmetries(1, 0) = -1;
  H_symmetries(1, 1) = 0;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandSpinDmn, typename KDmn>
void square_lattice<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KDmn>>& H_0) {
  static_assert(util::Length<typename BandSpinDmn::this_type>::value == 2,
                "BandSpinDmn has two subdomains, the band domain and the spin domain.");

  using BandDmn = typename util::TypeAt<0, typename BandSpinDmn::template domain_typelist<0>>::type;
  using SpinDmn = typename util::TypeAt<0, typename BandSpinDmn::template domain_typelist<1>>::type;

  if (BandDmn::get_size() != BANDS)
    throw std::logic_error("Square lattice has one band.");
  if (SpinDmn::get_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t = parameters.get_t();
  const auto t_prime = parameters.get_t_prime();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const auto val =
        -2. * t * (std::cos(k[0]) + std::cos(k[1])) - 4. * t_prime * std::cos(k[0]) * std::cos(k[1]);

    H_0(0, 0, 0, 0, k_ind) = val;
    H_0(0, 1, 0, 1, k_ind) = val;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
