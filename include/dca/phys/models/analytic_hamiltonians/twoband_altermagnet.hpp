// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas A. Maier (maierta@ornl.gov)
//
// Implementation of the Satoshi Okamoto's 2-orbital altermagnet Hubbard model

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TWOBAND_ALTERMAGNET_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TWOBAND_ALTERMAGNET_HPP

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

// TODO: the symmetry of this system must be checked.
template <typename SymmetryGroup>
class Altermagnet {
public:
  static constexpr bool complex_g0 = false;
  static constexpr bool spin_symmetric = true;

  using LDA_point_group = domains::no_symmetry<2>;
  using DCA_point_group = SymmetryGroup;
  // typedef PointGroupType DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 2;

  static double* initializeRDCABasis();
  // static double* initializeKDCABasis();

  static double* initializeRLDABasis();
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
int Altermagnet<PointGroupType>::transformationSignOfR(int b1 [[maybe_unused]],
                                                         int b2 [[maybe_unused]],
                                                         int s [[maybe_unused]]) {
  return 1;  // TODO: FIXME
}

template <typename PointGroupType>
int Altermagnet<PointGroupType>::transformationSignOfK(int b1 [[maybe_unused]],
                                                         int b2 [[maybe_unused]],
                                                         int s [[maybe_unused]]) {
  return 1;  // TODO: FIXME
}

template <typename PointGroupType>
double* Altermagnet<PointGroupType>::initializeRDCABasis() {
  static std::array<double, 4> basis{1.0, 1.0, 1.0, -1.0};
  return basis.data();
}

template <typename PointGroupType>
double* Altermagnet<PointGroupType>::initializeRLDABasis() {
  static std::array<double, 4> basis{1.0, 1.0, 1.0, -1.0};
  return basis.data();
}

template <typename PointGroupType>
std::vector<int> Altermagnet<PointGroupType>::flavors() {
  return {0, 1};
}

template <typename PointGroupType>
std::vector<std::vector<double>> Altermagnet<PointGroupType>::aVectors() {
  return {{0, 0}, {1.0, 0}};
}

template <typename PointGroupType>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> Altermagnet<
    PointGroupType>::orbitalPermutations() {
  return {};
}

template <typename PointGroupType>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void Altermagnet<PointGroupType>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Altermagnet Hubbard model has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // on-site Hubbard interaction

  H_interaction = 0.;

  for (int i = 0; i < BANDS; i++) {
    H_interaction(i, 0, i, 1, origin) = U;
    H_interaction(i, 1, i, 0, origin) = U;
  }
}

template <typename PointGroupType>
template <class domain>
void Altermagnet<PointGroupType>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

}

template <typename PointGroupType>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void Altermagnet<PointGroupType>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Altermagnet Hubbard model has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t0 = parameters.get_t0();
  const auto t1 = parameters.get_t1();
  const auto t2 = parameters.get_t2();

  H_0 = ScalarType(0);

  const ScalarType I(0, 1.);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    for (int s = 0; s < 2; s++) {
      H_0(0, s, 0, s, k_ind) = -2. * t1 * std::cos(2.0 * k[0]) - 2. * t2 * std::cos(2.0 * k[1]);
      H_0(1, s, 1, s, k_ind) = -2. * t2 * std::cos(2.0 * k[0]) - 2. * t1 * std::cos(2.0 * k[1]);
      H_0(0, s, 1, s, k_ind) = -2. * t0 * (std::cos(k[0]) + std::cos(k[1]));
      H_0(1, s, 0, s, k_ind) = -2. * t0 * (std::cos(k[0]) + std::cos(k[1]));
    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TWOBAND_ALTERMAGNET_HPP
