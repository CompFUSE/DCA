// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides a function to set a random vertex in the Monte Carlo solver.
//
// TODO: - Rename class.
//       - Const correctness.

#ifndef DCA_PHYS_MODELS_GENERAL_INTERACTION_HPP
#define DCA_PHYS_MODELS_GENERAL_INTERACTION_HPP

#include <array>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/convert.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename parameters_type>
class general_interaction {
public:
  // Initializes a (CT-AUX) vertex pair with a random pair of correlated spin-orbitals i.e. two
  // spin-orbitals with non-vanishing interaction.
  template <typename vertex_pair_type, typename rng_type, typename ScalarType, typename BandDmn,
            typename SpinDmn, typename RDmn>
  static void set_vertex(
      vertex_pair_type& vertex, const parameters_type& parameters, rng_type& rng,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>&
          H_interaction);

  // Creates and returns a vector of the linear indices of the correlated orbitals, i.e. pairs of
  // spin-orbitals with non-vanishing interaction. The linear index is with respect to the product
  // domain <band_1, spin_1, band_2, spin_2, site_1, site_2>.
  // Only orbitals set as interacting in the parameters are considered.
  template <typename ScalarType, typename BandDmn, typename SpinDmn, typename RDmn>
  static std::vector<int> make_correlated_orbitals(
      const parameters_type& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>&
          H_interaction);

private:
  // Creates and returns an object (singleton) of type spin_orbital_pair_domain, which is a product
  // domain of two band domains, two spin domains and two (real space) site domains.
  // The order of the subdomains is chosen to be: band_1, spin_1, band_2, spin_2, site_1, site_2.
  template <typename BandDmn, typename SpinDmn, typename RDmn>
  static func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, RDmn, RDmn>& get_spin_orbital_pair_domain() {
    static func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, RDmn, RDmn> spin_orbital_pair_domain;
    return spin_orbital_pair_domain;
  }
};

template <typename parameters_type>
template <typename vertex_pair_type, typename rng_type, typename ScalarType, typename BandDmn,
          typename SpinDmn, typename RDmn>
void general_interaction<parameters_type>::set_vertex(
    vertex_pair_type& vertex, const parameters_type& parameters, rng_type& rng,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>&
        H_interaction) {
  // Create the vector of correlated spin-orbitals once.
  static const std::vector<int> correlated_orbitals =
      general_interaction<parameters_type>::make_correlated_orbitals(parameters, H_interaction);

  // Get a random pair of correlated spin-orbitals.
  const int pos = rng() * correlated_orbitals.size();
  // This is instead of crashing on segv when module is unknown.
  // \todo catch earlier in release as well.
  assert( pos < correlated_orbitals.size() && "It is likely you have specified an unknown model" );
  const int lin_ind = correlated_orbitals[pos];

  std::array<int, 6> sub_ind;  // [0]=b1, [1]=s1, [2]=b2, [3]=s2, [4]=r1, [5]=r2
  get_spin_orbital_pair_domain<BandDmn, SpinDmn, RDmn>().linind_2_subind(lin_ind, sub_ind.data());

  // Set the vertex properties.
  vertex.get_bands().first = sub_ind[0];
  vertex.get_bands().second = sub_ind[2];

  vertex.get_e_spins().first = SpinDmn::get_elements()[sub_ind[1]];
  vertex.get_e_spins().second = SpinDmn::get_elements()[sub_ind[3]];

  vertex.get_spin_orbitals().first =
      domains::convert<int, func::dmn_variadic<BandDmn, SpinDmn>>::spin_orbital(
          vertex.get_bands().first, vertex.get_e_spins().first);
  vertex.get_spin_orbitals().second =
      domains::convert<int, func::dmn_variadic<BandDmn, SpinDmn>>::spin_orbital(
          vertex.get_bands().second, vertex.get_e_spins().second);

  vertex.get_r_sites().first = sub_ind[4];
  vertex.get_r_sites().second = sub_ind[5];
}

template <typename parameters_type>
template <typename ScalarType, typename BandDmn, typename SpinDmn, typename RDmn>
std::vector<int> general_interaction<parameters_type>::make_correlated_orbitals(
    const parameters_type& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>&
        H_interaction) {
  std::vector<int> correlated_orbitals;

  const std::vector<int>& interacting_orbitals = parameters.get_interacting_orbitals();

  for (int r_j = 0; r_j < RDmn::dmn_size(); ++r_j) {
    for (int r_i = 0; r_i < RDmn::dmn_size(); ++r_i) {
      const int delta_r = RDmn::parameter_type::subtract(r_j, r_i);  // delta_r = r_i - r_j

      for (int s_j = 0; s_j < SpinDmn::dmn_size(); ++s_j) {
        for (auto b_j : interacting_orbitals) {
          for (int s_i = 0; s_i < SpinDmn::dmn_size(); ++s_i) {
            for (auto b_i : interacting_orbitals) {
              // This 1.e-3 seems like a very magic number!
              if (std::abs(H_interaction(b_i, s_i, b_j, s_j, delta_r)) > 1.e-3) {
                int linear_index = get_spin_orbital_pair_domain<BandDmn, SpinDmn, RDmn>()(
                    b_i, s_i, b_j, s_j, r_i, r_j);
                correlated_orbitals.push_back(linear_index);
              }
            }
          }
        }
      }
    }
  }

  return correlated_orbitals;
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_GENERAL_INTERACTION_HPP
