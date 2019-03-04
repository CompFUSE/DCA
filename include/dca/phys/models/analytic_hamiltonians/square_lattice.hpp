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

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP

#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/util.hpp"

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

  // Initializes the interaction part of the real space Hubbard Hamiltonian.
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

  // Initializes the hopping matrix for the case of broken translational invariance in the cluster.
  // Then, the hopping matrix is a function of real space cluster vectors I and J, and a real space
  // superlattice vector d_tilde (= i_tilde - j_tilde, translational invariance of the
  // superlattice).
  template <typename RLatticeDmn, typename ParametersType, typename ScalarType,
            typename RClusterDmn, typename RSuperlatticeDmn>
  static void initializeNonTranslationalInvariantHoppingMatrix(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<RClusterDmn, RClusterDmn, RSuperlatticeDmn>>&
          t_IJ_d_tilde);
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

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void square_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
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

  util::initializeSingleBandHint(parameters, nn_vec, H_interaction);
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
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void square_lattice<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
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

template <typename point_group_type>
template <typename RLatticeDmn, typename ParametersType, typename ScalarType, typename RClusterDmn,
          typename RSuperlatticeDmn>
void square_lattice<point_group_type>::initializeNonTranslationalInvariantHoppingMatrix(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<RClusterDmn, RClusterDmn, RSuperlatticeDmn>>&
        t_IJ_d_tilde) {
  if (parameters.get_t_prime() != 0.)
    throw std::logic_error("Next nearest neighbor hopping has not been implemented yet.");

  const auto& basis = RClusterDmn::parameter_type::get_basis_vectors();
  // The four directions of nearest neighbor hopping.
  const auto a0 = basis[0];
  const auto min_a0 = math::util::scale(-1., a0);
  const auto a1 = basis[1];
  const auto min_a1 = math::util::scale(-1., a1);

  const auto& lattice_superbasis = RLatticeDmn::parameter_type::get_super_basis_vectors();
  // To apply the periodic boundary conditions of the lattice.
  const auto b0 = lattice_superbasis[0];
  const auto b1 = lattice_superbasis[1];
  const auto b0_plus_b1 = math::util::add(b0, b1);

  const auto& cluster = RClusterDmn::get_elements();
  const auto& superlattice = RSuperlatticeDmn::get_elements();

  const auto t = parameters.get_t();

  for (std::size_t I = 0; I < cluster.size(); ++I) {
    const auto& I_vec = cluster[I];
    for (std::size_t J = 0; J < cluster.size(); ++J) {
      const auto& J_vec = cluster[J];
      for (std::size_t d_tilde = 0; d_tilde < superlattice.size(); ++d_tilde) {
        const auto& d_tilde_vec = superlattice[d_tilde];

        // i = I + d_tilde, j = J
        const auto i_min_j = math::util::subtract(J_vec, math::util::add(I_vec, d_tilde_vec));

        // Check whether i and j are nearest neighbors.
        if (math::util::isSameVector(i_min_j, a0) || math::util::isSameVector(i_min_j, min_a0) ||
            math::util::isSameVector(i_min_j, a1) || math::util::isSameVector(i_min_j, min_a1))
          t_IJ_d_tilde(I, J, d_tilde) = t;

        // Take into account the periodic boundary conditions of the lattice.
        // Translate i-j by b0.
        const auto i_min_j_min_b0 = math::util::subtract(b0, i_min_j);
        if (math::util::isSameVector(i_min_j_min_b0, a0) ||
            math::util::isSameVector(i_min_j_min_b0, min_a0) ||
            math::util::isSameVector(i_min_j_min_b0, a1) ||
            math::util::isSameVector(i_min_j_min_b0, min_a1))
          t_IJ_d_tilde(I, J, d_tilde) = t;

        // Translate i-j by b1.
        const auto i_min_j_min_b1 = math::util::subtract(b1, i_min_j);
        if (math::util::isSameVector(i_min_j_min_b1, a0) ||
            math::util::isSameVector(i_min_j_min_b1, min_a0) ||
            math::util::isSameVector(i_min_j_min_b1, a1) ||
            math::util::isSameVector(i_min_j_min_b1, min_a1))
          t_IJ_d_tilde(I, J, d_tilde) = t;

        // Translate i-j by b0 + b1.
        const auto i_min_j_min_b0_min_b1 = math::util::subtract(b0_plus_b1, i_min_j);
        if (math::util::isSameVector(i_min_j_min_b0_min_b1, a0) ||
            math::util::isSameVector(i_min_j_min_b0_min_b1, min_a0) ||
            math::util::isSameVector(i_min_j_min_b0_min_b1, a1) ||
            math::util::isSameVector(i_min_j_min_b0_min_b1, min_a1))
          t_IJ_d_tilde(I, J, d_tilde) = t;
      }
    }
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
