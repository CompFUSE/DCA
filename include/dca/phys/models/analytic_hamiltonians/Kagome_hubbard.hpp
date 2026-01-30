// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas A. Maier (maierta@ornl.gov)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of the Kagome Hubbard model

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KAGOME_HUBBARD_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KAGOME_HUBBARD_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

// TODO: the symmetry of this system must be checked.
template <typename SymmetryGroup>
class KagomeHubbard {
public:
  static constexpr bool complex_g0 = false;
  static constexpr bool spin_symmetric = true;

  using LDA_point_group = domains::no_symmetry<2>;
  using DCA_point_group = SymmetryGroup;
  // typedef PointGroupType DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 3;

  const static double* initializeRDCABasis();
  // static double* initializeKDCABasis();

  const static double* initializeRLDABasis();
  // static double* initializeKLDABasis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Rotations of pi/2 are an anti-symmetry on the band off-diagonal.
  static int transformationSignOfR(int b1, int b2, int s);
  static int transformationSignOfK(int b1, int b2, int s);

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
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
};

template <typename PointGroupType>
int KagomeHubbard<PointGroupType>::transformationSignOfR(int b1 [[maybe_unused]],
                                                         int b2 [[maybe_unused]],
                                                         int s [[maybe_unused]]) {
  return 1;  // TODO: FIXME
}

template <typename PointGroupType>
int KagomeHubbard<PointGroupType>::transformationSignOfK(int b1 [[maybe_unused]],
                                                         int b2 [[maybe_unused]],
                                                         int s [[maybe_unused]]) {
  return 1;  // TODO: FIXME
}

template <typename PointGroupType>
const double* KagomeHubbard<PointGroupType>::initializeRDCABasis() {
  static std::array<double, 4> basis{std::cos(M_PI / 3.), std::sin(M_PI / 3.), std::cos(M_PI / 3.),
                                     std::sin(-M_PI / 3.)};
  // static std::array<double, 4> basis{1, 0, 0.5, std::sin(M_PI / 3.)};
  return basis.data();
}
// template <typename PointGroupType>
// double* KagomeHubbard<PointGroupType>::initializeKDCABasis() {
//   static std::array<double, 4> basis{2 * M_PI, -M_PI/std::sin(M_PI/3.), 0, 2. *
//   M_PI/std::sin(M_PI/3.)}; return basis.data();
// }

template <typename PointGroupType>
const double* KagomeHubbard<PointGroupType>::initializeRLDABasis() {
  static std::array<double, 4> basis{std::cos(M_PI / 3.), std::sin(M_PI / 3.), std::cos(M_PI / 3.),
                                     std::sin(-M_PI / 3.)};
  // static std::array<double, 4> basis{1, 0, 0, std::sin(M_PI / 3.)};
  return basis.data();
}

// template <typename PointGroupType>
// double* KagomeHubbard<PointGroupType>::initializeKLDABasis() {
//   static std::array<double, 4> basis{2 * M_PI, -M_PI/std::sin(M_PI/3.), 0, 2. *
//   M_PI/std::sin(M_PI/3.)}; return basis.data();
// }

template <typename PointGroupType>
std::vector<int> KagomeHubbard<PointGroupType>::flavors() {
  return {0, 1, 2};
}

template <typename PointGroupType>
std::vector<std::vector<double>> KagomeHubbard<PointGroupType>::aVectors() {
  return {{0, 0}, {0.5, 0}, {0.25, std::sin(M_PI / 3.) / 2.}};
}

template <typename PointGroupType>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> KagomeHubbard<
    PointGroupType>::orbitalPermutations() {
  return {};
}

template <typename PointGroupType>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void KagomeHubbard<PointGroupType>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Kagome Hubbard model has three bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // on-site Hubbard interaction
  const double V = parameters.get_V();  // Nearest-neighbor interaction

  // There are 6 different nearest neighbor (nn) pairs along a1, a2, a1+a2 and their negative vectors
  const std::vector<typename RDmn::parameter_type::element_type>& basis =
      RDmn::parameter_type::get_basis_vectors();

  assert(basis.size() == 2);

  // Pre-compute basis vectors and orbital pairs
  std::vector<typename RDmn::parameter_type::element_type> nn_vec;
  nn_vec.reserve(3);
  nn_vec.push_back(basis[0]);
  nn_vec.push_back(basis[1]);
  nn_vec.emplace_back(basis[0]);
  nn_vec[2][0] += basis[1][0];
  nn_vec[2][1] += basis[1][1];

  // Orbital pairs corresponding to each nearest neighbor direction
  constexpr std::array<std::array<int, 2>, 3> orb_pairs = {{
    {{2, 0}}, {{1, 2}}, {{1, 0}}
  }};

  const std::vector<typename RDmn::parameter_type::element_type>& super_basis =
      RDmn::parameter_type::get_super_basis_vectors();
  const std::vector<typename RDmn::parameter_type::element_type>& elements =
      RDmn::parameter_type::get_elements();

  std::vector<int> nn_index;
  std::vector<std::vector<int>> nn_orbs;
  nn_index.reserve(6);
  nn_orbs.reserve(6);

  // Process all nearest neighbor vectors in a loop
  for (int vec_idx = 0; vec_idx < 3; ++vec_idx) {
    const auto nn_vec_translated =
        domains::cluster_operations::translate_inside_cluster(nn_vec[vec_idx], super_basis);
    const int r = domains::cluster_operations::index(nn_vec_translated, elements, domains::BRILLOUIN_ZONE);

    if (r != origin) {
      const int minus_r = RDmn::parameter_type::subtract(r, origin);

      nn_index.push_back(r);
      nn_orbs.push_back({orb_pairs[vec_idx][0], orb_pairs[vec_idx][1]});

      // Add negative direction if distinct
      if (r != minus_r) {
        nn_index.push_back(minus_r);
        nn_orbs.push_back({orb_pairs[vec_idx][1], orb_pairs[vec_idx][0]});
      }
    }
  }

  H_interaction = 0.;

// Onsite U interaction
  for (int i = 0; i < BANDS; i++) {
    H_interaction(i, 0, i, 1, origin) = U;
    H_interaction(i, 1, i, 0, origin) = U;
  }

  // Nearest-neighbor V interaction inside unit cell
  constexpr std::array<std::array<int, 2>, 3> intra_cell_pairs = {{
    {{0, 1}}, {{1, 2}}, {{0, 2}}
  }};

  for (const auto& pair : intra_cell_pairs) {
    for (int s = 0; s < 2; ++s) {
      H_interaction(pair[0], s, pair[1], s, origin) = V;
      H_interaction(pair[0], s, pair[1], 1-s, origin) = V;
      H_interaction(pair[1], s, pair[0], s, origin) = V;
      H_interaction(pair[1], s, pair[0], 1-s, origin) = V;
    }
  }

  // Nearest-neighbor V interaction outside unit cell
  const auto nn_size = nn_index.size();
  for (std::size_t i = 0; i < nn_size; ++i) {
    for (int s = 0; s < 2; ++s) {
      H_interaction(nn_orbs[i][0], s, nn_orbs[i][1], s, nn_index[i]) = V;
      H_interaction(nn_orbs[i][0], s, nn_orbs[i][1], 1-s, nn_index[i]) = V;
    }
  }

}




template <typename PointGroupType>
template <class domain>
void KagomeHubbard<PointGroupType>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  //  H_symmetry(i, s1, j, s2)
  //  H_symmetries(0, 0, 0, 0) = 0; // at b0, G of spin 0 or 1 has the same values.
  //  H_symmetries(0, 1, 0, 1) = 0;

  //  H_symmetries(1, 0, 1, 0) = 1;
  //  H_symmetries(1, 1, 1, 1) = 1; // at i, G of spin 0 or 1 has the same values.
}

template <typename PointGroupType>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void KagomeHubbard<PointGroupType>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Kagome Hubbard model has three bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t = parameters.get_t();

  H_0 = ScalarType(0);

  const ScalarType I(0, 1.);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    for (int s = 0; s < 2; s++) {
      H_0(0, s, 1, s, k_ind) = -2. * t * std::cos(0.5 * k[0]);
      H_0(0, s, 2, s, k_ind) = -2. * t * std::cos(0.25 * k[0] + 0.25 * sqrt(3.) * k[1]);
      H_0(1, s, 2, s, k_ind) = -2. * t * std::cos(0.25 * k[0] - 0.25 * sqrt(3.) * k[1]);

      H_0(1, s, 0, s, k_ind) = -2. * t * std::cos(0.5 * k[0]);
      H_0(2, s, 0, s, k_ind) = -2. * t * std::cos(0.25 * k[0] + 0.25 * sqrt(3.) * k[1]);
      H_0(2, s, 1, s, k_ind) = -2. * t * std::cos(0.25 * k[0] - 0.25 * sqrt(3.) * k[1]);
    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_KAGOME_HUBBARD_HPP
