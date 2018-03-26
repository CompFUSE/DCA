// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename PointGroupType>
class FeAsLattice : public bilayer_lattice<PointGroupType> {
public:
  using BaseClass = bilayer_lattice<PointGroupType>;
  constexpr static int BANDS = BaseClass::BANDS;
  constexpr static int DIMENSION = BaseClass::DIMENSION;

  // Initialize the non interacting Hamiltonian.
  template <typename Parameters, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initialize_H_0(
      const Parameters& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  // Initializes the interaction Hamiltonian  density-density local term.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
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
void FeAsLattice<PointGroupType>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();
  const auto t = parameters.get_t();

  H_0 = ScalarType(0);
  using std::cos;
  using std::sin;

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const double kx(k_vecs[k_ind][0]), ky(k_vecs[k_ind][1]);
    // Same band interaction.
    const double eps_0 = -2 * t[0] * cos(kx) - 2 * t[1] * cos(ky) - 4 * t[2] * cos(kx) * cos(ky);
    const double eps_1 = -2 * t[1] * cos(kx) - 2 * t[0] * cos(ky) - 4 * t[2] * cos(kx) * cos(ky);
    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_ind) = eps_0;
      H_0(1, s, 1, s, k_ind) = eps_1;
    }
    // Off band interaction
    const double eps_01 = -4 * sin(kx) * sin(ky);
    for (int s = 0; s < 2; ++s)
      H_0(0, s, 1, s, k_ind) =  H_0(1, s, 0, s, k_ind) = eps_01;
  }
}

template <typename PointGroupType>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void FeAsLattice<PointGroupType>::initialize_H_interaction(
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
          // Coulomb repulsion and contribution from S*S interaction.
          if (b1 == b2 and s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;
          else if (b1 != b2 and s1 == s2)
            H_interaction(b1, s1, b2, s2, origin) = V + J;
          else if (b1 != b2 and s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = V - J;
        }
}

template <typename PointGroupType>
template <typename Nu, typename RDmn, typename parameters_type>
void FeAsLattice<PointGroupType>::initializeNonDensityInteraction(
    func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, RDmn>>& non_density_interaction,
    const parameters_type& parameters) {
  const double J = parameters.get_J();
  const Nu nu;  // band-spin domain.
  constexpr int up(0), down(1);

  non_density_interaction = 0.;
  for (int b1 = 0; b1 < BANDS; b1++)
    for (int b2 = 0; b2 < BANDS; b2++) {
      if (b1 == b2)
        continue;
      // spin-flip interaction equivalent to  S^+S^- and S^-S^+ interactions.
      non_density_interaction(nu(b1, up), nu(b2, up), nu(b2, down), nu(b1, down), 0) = -J;
    }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_LATTICE_HPP
