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
// Implementation of  the two-orbital model described in "Two pairing domes as Cu2+ varies to Cu3+."
// https://link.aps.org/doi/10.1103/PhysRevB.99.224515

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TWOBAND_CU_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TWOBAND_CU_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename PointGroup>
class TwoBandCu {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef PointGroup DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 2;

  // Rotations of pi/2 are an anti-symmetry on the band off-diagonal.
  static int transformationSign(int b1, int b2, int s);

  static int transformationSignOfR(int b1, int b2, int s) {
    return transformationSign(b1, b2, s);
  }
  static int transformationSignOfK(int b1, int b2, int s) {
    return transformationSign(b1, b2, s);
  }

  static double* initializeRDCABasis();

  static double* initializeKDCABasis();

  static double* initializeRLDABasis();

  static double* initializeKLDABasis();

  static std::vector<int> flavors();

  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class HType, class Parameters>
  static void initializeNonDensityInteraction(HType& interaction, const Parameters& pars);

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

template <typename PointGroup>
int TwoBandCu<PointGroup>::transformationSign(int b1, int b2, int s) {
  if (!std::is_same<PointGroup, domains::D4>::value)
    return 1;

  if (b1 == b2)
    return 1;

  constexpr bool symmetrize_off_diagonal = true;

  if (symmetrize_off_diagonal) {
    // TODO: generalize.

    const bool is_odd_rotation = s == 1 || s == 3;
    const bool is_odd_reflection = s == 4 || s == 6;

    return is_odd_reflection || is_odd_rotation ? -1 : 1;
  }

  else {
    // TODO: generalize.
    const bool is_identity = s == 0;
    return is_identity ? 1 : 0;
  }
}

template <typename PointGroup>
double* TwoBandCu<PointGroup>::initializeRDCABasis() {
  static std::array<double, 4> r_DCA{1, 0, 0, 1};
  return r_DCA.data();
}

template <typename PointGroup>
double* TwoBandCu<PointGroup>::initializeKDCABasis() {
  static std::array<double, 4> k_DCA{2 * M_PI, 0, 0, 2 * M_PI};
  return k_DCA.data();
}

template <typename PointGroup>
double* TwoBandCu<PointGroup>::initializeRLDABasis() {
  static std::array<double, 4> basis{1, 0, 0, 1};
  return basis.data();
}

template <typename PointGroup>
double* TwoBandCu<PointGroup>::initializeKLDABasis() {
  static std::array<double, 4> basis{2 * M_PI, 0, 0, 2 * M_PI};
  return basis.data();
}

template <typename PointGroup>
std::vector<int> TwoBandCu<PointGroup>::flavors() {
  return {0, 1};
}

template <typename PointGroup>
std::vector<std::vector<double>> TwoBandCu<PointGroup>::aVectors() {
  return {{0, 0}, {0, 0}};
}

template <typename PointGroup>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> TwoBandCu<PointGroup>::orbitalPermutations() {
  return {};
}

template <typename PointGroup>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void TwoBandCu<PointGroup>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("2 orbital lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // Intra-orbital interaction.
  const double V = parameters.get_V();  // Inter-orbital interaction.
  const double J = parameters.get_J();  // Spin-spin term.
  // const double V_prime = parameters.get_V_prime();  // Different orbital, same spin.

  H_interaction = 0.;

  // Note: the solvers expect a double counted Hamiltonian.
  for (int b1 = 0; b1 < BANDS; b1++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int b2 = 0; b2 < BANDS; b2++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (b1 == b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;

          if (b1 != b2)
            H_interaction(b1, s1, b2, s2, origin) = V - (s1 == s2) * J;
        }
      }
    }
  }
}

template <typename PointGroup>
template <class HType, class Parameters>
void TwoBandCu<PointGroup>::initializeNonDensityInteraction(HType& interaction,
                                                            const Parameters& pars) {
  const double J = pars.get_J();
  const double Jp = pars.get_Jp();

  constexpr int up(0), down(1);

  interaction = 0.;
  for (int b1 = 0; b1 < BANDS; b1++)
    for (int b2 = 0; b2 < BANDS; b2++) {
      if (b1 == b2)
        continue;
      interaction(b1, up, b2, up, b2, down, b1, down, 0) = J;
      interaction(b1, up, b2, up, b1, down, b2, down, 0) = Jp;
    }
}
template <typename PointGroup>
template <class domain>
void TwoBandCu<PointGroup>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
}

template <typename PointGroup>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void TwoBandCu<PointGroup>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("TwoBandCu lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t00 = parameters.get_t00();
  const auto tp00 = parameters.get_tp00();
  const auto tpp00 = parameters.get_tpp00();

  const auto t11 = parameters.get_t11();
  const auto tp11 = parameters.get_tp11();
  const auto tpp11 = parameters.get_tpp11();

  const auto t01 = parameters.get_t01();
  // const auto tp01   =  parameters.get_tp01() ;
  const auto tpp01 = parameters.get_tpp01();

  const auto delta_E = parameters.get_delta_E();

  H_0 = 0;

  for (int k_id = 0; k_id < k_vecs.size(); ++k_id) {
    const auto& k = k_vecs[k_id];

    // const auto val00 = 2.*t00*(std::cos(k[0])+std::cos(k[1]));
    const auto val00 = 2. * t00 * (std::cos(k[0]) + std::cos(k[1])) +
                       4. * tp00 * std::cos(k[0]) * std::cos(k[1]) +
                       2. * tpp00 * (std::cos(2. * k[0]) + std::cos(2. * k[1]));

    // const auto val11 = 2.*t11*(std::cos(k[0])+std::cos(k[1]));
    const auto val11 = delta_E + 2. * t11 * (std::cos(k[0]) + std::cos(k[1])) +
                       4. * tp11 * std::cos(k[0]) * std::cos(k[1]) +
                       2. * tpp11 * (std::cos(2. * k[0]) + std::cos(2. * k[1]));

    // const auto val01 = t01;
    const auto val01 = 2. * t01 * (std::cos(k[0]) - std::cos(k[1])) +
                       2. * tpp01 * (std::cos(2. * k[0]) - std::cos(2. * k[1]));

    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_id) = val00;
      H_0(1, s, 1, s, k_id) = val11;

      H_0(0, s, 1, s, k_id) = val01;
      H_0(1, s, 0, s, k_id) = val01;
    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_LATTICE_HPP
