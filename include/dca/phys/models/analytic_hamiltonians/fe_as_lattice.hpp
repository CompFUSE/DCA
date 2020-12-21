// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Bilayer lattice with spin-spin interaction. See "model_parameters_fe_as.hpp".

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_LATTICE_HPP

#include <cmath>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

struct FeAsPointGroup {
  using point_group_type_list =
      dca::util::Typelist<domains::identity_group_operation<2>, domains::Cn_2D<2, 4>>;
};

template <typename /*PointGroupType*/>
class FeAsLattice : public bilayer_lattice<FeAsPointGroup> {
public:
  using BaseClass = bilayer_lattice<FeAsPointGroup>;
  constexpr static int BANDS = BaseClass::BANDS;
  constexpr static int DIMENSION = BaseClass::DIMENSION;

  // Initialize the non interacting Hamiltonian.
  template <typename Parameters, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const Parameters& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  // Initializes the interaction Hamiltonian  density-density local term.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  // Initializes the interaction Hamiltonian non density-density local term.
  template <typename Nu, typename RDmn, typename parameters_type>
  static void initializeNonDensityInteraction(
      func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
      const parameters_type& parameters);
};

template <typename PointGroupType>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void FeAsLattice<PointGroupType>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto t1 = parameters.get_t1();
  const auto t2 = parameters.get_t2();
  const auto t3 = parameters.get_t3();
  const auto t4 = parameters.get_t4();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = KDmn::get_elements()[k_ind];

    const double cx = std::cos(k[0]);
    const double cy = std::cos(k[1]);

    const double val_00 = -2 * t1 * cx - 2 * t2 * cy - 4 * t3 * cx * cy;
    const double val_11 = -2 * t2 * cx - 2 * t1 * cy - 4 * t3 * cx * cy;
    const double val_01 = -4 * t4 * std::sin(k[0]) * std::sin(k[1]);

    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_ind) = val_00;
      H_0(1, s, 1, s, k_ind) = val_11;

      H_0(0, s, 1, s, k_ind) = val_01;
      H_0(1, s, 0, s, k_ind) = val_01;
    }
  }
}

template <typename PointGroupType>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void FeAsLattice<PointGroupType>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // Same band interaction.
  const double V = parameters.get_V();  // Different band interaction.
  const double J = parameters.get_J();  // Spin-spin interaction.
  H_interaction = 0.;
  for (int b1 = 0; b1 < BANDS; ++b1)
    for (int b2 = 0; b2 < BANDS; ++b2)
      for (int s1 = 0; s1 < 2; ++s1)
        for (int s2 = 0; s2 < 2; s2++) {
          // Coulomb repulsion and contribution from -J S_z*S_z interaction.
          if (b1 == b2 and s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;
          else if (b1 != b2 and s1 == s2)
            H_interaction(b1, s1, b2, s2, origin) = V - J;
          else if (b1 != b2 and s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = V + J;
        }
}

template <typename PointGroupType>
template <typename Nu, typename RDmn, typename parameters_type>
void FeAsLattice<PointGroupType>::initializeNonDensityInteraction(
    func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
    const parameters_type& parameters) {
  const double J = parameters.get_J();
  const double Jp = parameters.get_Jp();
  const Nu nu;  // band-spin domain.
  constexpr int up(0), down(1);

  non_density_interaction = 0.;
  for (int b1 = 0; b1 < BANDS; b1++)
    for (int b2 = 0; b2 < BANDS; b2++) {
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

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_LATTICE_HPP
