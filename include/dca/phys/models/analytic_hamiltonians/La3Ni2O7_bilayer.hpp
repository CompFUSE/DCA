// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas Maier maierta@ornl.gov
//
// Two-orbital bilayer model for La3Ni2O7

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_La3Ni2O7_BILAYER_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_La3Ni2O7_BILAYER_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"
#include "dca/phys/models/traits.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class La3Ni2O7_bilayer {
public:
  static constexpr bool complex_g0 = false;
  static constexpr bool spin_symmetric = true;

  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 4;

  constexpr static int transformationSignOfR(int, int, int) {
    return 1;
  }
  constexpr static int transformationSignOfK(int, int, int) {
    return 1;
  }

  static const double* initializeRDCABasis();
  static const double* initializeKDCABasis();

  static const double* initializeRLDABasis();
  static const double* initializeKLDABasis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  // Initializes the interaction Hamiltonian non density-density local term.
  template <typename Scalar, typename Parameters>
  static void initializeNonDensityInteraction(
					      NonDensityIntHamiltonian<Scalar, Parameters>& non_density_interaction, const Parameters& parameters);

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

template <typename point_group_type>
const double* La3Ni2O7_bilayer<point_group_type>::initializeRDCABasis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
const double* La3Ni2O7_bilayer<point_group_type>::initializeKDCABasis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
const double* La3Ni2O7_bilayer<point_group_type>::initializeRLDABasis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
const double* La3Ni2O7_bilayer<point_group_type>::initializeKLDABasis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> La3Ni2O7_bilayer<point_group_type>::flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;
  flavors[2] = 2;
  flavors[3] = 3;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> La3Ni2O7_bilayer<point_group_type>::aVectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> La3Ni2O7_bilayer<
    point_group_type>::orbitalPermutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void La3Ni2O7_bilayer<point_group_type>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has four bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  H_interaction = 0.;

  const double U = parameters.get_U();  // Intra-orbital interaction.
  const double V = parameters.get_V();  // Inter-orbital interaction.
  const double J = parameters.get_J();  // Spin-spin interaction.

  // Basis: 0=dx2-y2,l1; 1=d3z2,l1; 2=dx2-y2,l1; 3=dz2,l2
  for (int b = 0; b < BANDS; ++b) {
        H_interaction(b, 0, b, 1, origin) = U;
        H_interaction(b, 1, b, 0, origin) = U;
  }
  for (int s = 0; s < 2; ++s) {
        H_interaction(0, s, 1, s, origin) = V - J;
        H_interaction(1, s, 0, s, origin) = V - J;
        H_interaction(2, s, 3, s, origin) = V - J;
        H_interaction(3, s, 2, s, origin) = V - J;

        H_interaction(0, s, 1, 1-s, origin) = V;
        H_interaction(1, s, 0, 1-s, origin) = V;
        H_interaction(2, s, 3, 1-s, origin) = V;
        H_interaction(3, s, 2, 1-s, origin) = V;
  }


}

template <typename point_group_type>
template <class domain>
void La3Ni2O7_bilayer<point_group_type>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  /* H_symmetries(0, 0, 0, 0) = 0; */
  /* H_symmetries(0, 1, 0, 1) = 0; */
  /**/
  /* H_symmetries(1, 0, 1, 0) = 1; */
  /* H_symmetries(1, 1, 1, 1) = 1; */
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void La3Ni2O7_bilayer<point_group_type>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t11 = parameters.get_t11();
  const auto t12 = parameters.get_t12();
  const auto t22 = parameters.get_t22();
  const auto t_perp_11 = parameters.get_t_perp_11();
  const auto t_perp_22 = parameters.get_t_perp_22();
  const auto Delta = parameters.get_Delta();

  H_0 = ScalarType(0);

  // Basis: 0=dx2-y2,l1; 1=d3z2,l1; 2=dx2-y2,l1; 3=dz2,l2
  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const auto val11 = -2. * t11 * (std::cos(k[0]) + std::cos(k[1]));
    const auto val22 = -2. * t22 * (std::cos(k[0]) + std::cos(k[1]));
    const auto val12 = 2. * t12 * (std::cos(k[0]) - std::cos(k[1]));

    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_ind) = val11 + Delta;
      H_0(1, s, 1, s, k_ind) = val22;
      H_0(2, s, 2, s, k_ind) = val11 + Delta;
      H_0(3, s, 3, s, k_ind) = val22;

      H_0(0, s, 1, s, k_ind) = val12;
      H_0(1, s, 0, s, k_ind) = val12;
      H_0(2, s, 3, s, k_ind) = val12;
      H_0(3, s, 2, s, k_ind) = val12;

      H_0(0, s, 2, s, k_ind) = -t_perp_11;
      H_0(2, s, 0, s, k_ind) = -t_perp_11;
      H_0(1, s, 3, s, k_ind) = -t_perp_22;
      H_0(3, s, 1, s, k_ind) = -t_perp_22;
    }
  }
}

template <typename PointGroupType>
template <typename Scalar, typename Parameters>
void La3Ni2O7_bilayer<PointGroupType>::initializeNonDensityInteraction(
								  NonDensityIntHamiltonian<Scalar, Parameters>& non_density_interaction, const Parameters& parameters) {
  const double J = parameters.get_J();
  const double Jp = parameters.get_J();
  const NuDmn nu;  // band-spin domain.
  constexpr int up(0), down(1);

  non_density_interaction = 0.;
  for (int b1 = 0; b1 < 2; b1++)
    for (int b2 = 0; b2 < 2; b2++) {
      if (b1 == b2)
        continue;
      // spin-flip interaction coming from the -J * S_b1^+S_b2^- Hamiltonian term.
      // Note: a factor of -1 comes from rearranging the fermion operators.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b2, down), nu(b1, down), 0) = J;
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b1, down), nu(b2, down), 0) = Jp;
    }
  for (int b1 = 2; b1 < 4; b1++)
    for (int b2 = 2; b2 < 4; b2++) {
      if (b1 == b2)
        continue;
      // spin-flip interaction coming from the -J * S_b1^+S_b2^- Hamiltonian term.
      // Note: a factor of -1 comes from rearranging the fermion operators.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b2, down), nu(b1, down), 0) = J;
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b1, down), nu(b2, down), 0) = Jp;
    }
}


}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_La3Ni2O7_BILAYER_HPP
