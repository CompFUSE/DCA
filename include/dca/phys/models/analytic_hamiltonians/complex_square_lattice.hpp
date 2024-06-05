// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Square lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_COMPLEX_SQUARE_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_COMPLEX_SQUARE_LATTICE_HPP

#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/util.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class complex_square_lattice {
public:
  static constexpr bool complex_g0 = true;
  static constexpr bool spin_symmetric = true;

  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static int DIMENSION = 2;
  const static int BANDS = 1;

  static const double* initializeRDCABasis();
  static const double* initializeRLDABasis();

  constexpr static int transformationSignOfR(int, int, int) {
    return 1;
  }
  constexpr static int transformationSignOfK(int, int, int) {
    return 1;
  }

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction part of the real space Hubbard Hamiltonian.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
				     func::function<typename parameters_type::Real, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initializeHSymmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are two-dimensional (access through index 0 and 1).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0WithQ(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
      func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0, typename KDmn::element_type& q);

};

template <typename point_group_type>
const double* complex_square_lattice<point_group_type>::initializeRDCABasis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.;
  r_DCA[1] = 0.;
  r_DCA[2] = 0.;
  r_DCA[3] = 1.;

  return r_DCA;
}

template <typename point_group_type>
const double* complex_square_lattice<point_group_type>::initializeRLDABasis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
std::vector<int> complex_square_lattice<point_group_type>::flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> complex_square_lattice<point_group_type>::aVectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> complex_square_lattice<
    point_group_type>::orbitalPermutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void complex_square_lattice<point_group_type>::initializeHInteraction(
								      func::function<typename parameters_type::Real, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Square lattice has one band.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const std::vector<typename RDmn::parameter_type::element_type>& basis =
      RDmn::parameter_type::get_basis_vectors();

  assert(basis.size() == 2);

  // There are two different nearest neighbor (nn) pairs: along the basis vector a1 and along the
  // basis vector a2.
  const std::vector<typename RDmn::parameter_type::element_type>& nn_vec(basis);

  // Compute indices of nearest neighbors (nn) w.r.t. origin.
  std::vector<int> nn_index;

  const std::vector<typename RDmn::parameter_type::element_type>& super_basis =
      RDmn::parameter_type::get_super_basis_vectors();
  const std::vector<typename RDmn::parameter_type::element_type>& elements =
      RDmn::parameter_type::get_elements();

  for (const auto& vec : nn_vec) {
    std::vector<double> nn_vec_translated =
        domains::cluster_operations::translate_inside_cluster(vec, super_basis);
    nn_index.push_back(
        domains::cluster_operations::index(nn_vec_translated, elements, domains::BRILLOUIN_ZONE));
  }

  H_interaction = 0;

  const double V = parameters.get_V();
  for (auto index : nn_index) {
    H_interaction(0, 0, 0, 1, index) = V;
    H_interaction(0, 1, 0, 0, index) = V;
  }

  // Nearest-neighbor same spin interaction
  const double V_prime = parameters.get_V_prime();
  for (auto index : nn_index) {
    H_interaction(0, 0, 0, 0, index) = V_prime;
    H_interaction(0, 1, 0, 1, index) = V_prime;
  }

  // Get the index of the origin (0,0).
  const int origin = RDmn::parameter_type::origin_index();

  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  const double U = parameters.get_U();
  H_interaction(0, 0, 0, 1, origin) = U;
  H_interaction(0, 1, 0, 0, origin) = U;

  //  util::initializeSingleBandHint(parameters, nn_vec, H_interactioneraction);
}

template <typename point_group_type>
template <class domain>
void complex_square_lattice<point_group_type>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries(0, 0) = 0;
  H_symmetries(0, 1) = -1;
  H_symmetries(1, 0) = -1;
  H_symmetries(1, 1) = 0;
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void complex_square_lattice<point_group_type>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  typename KDmn::element_type default_q;
  initializeH0WithQ(parameters, H_0, default_q);
}
  
template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void complex_square_lattice<point_group_type>::initializeH0WithQ(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0, typename KDmn::element_type& q [[maybe_unused]]) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Square lattice has one band.");
  if (SpinDmn::dmn_size() != 2)
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

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_COMPLEX_SQUARE_LATTICE_HPP
