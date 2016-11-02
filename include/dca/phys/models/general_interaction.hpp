// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides a function to set a vertex in the Monte Carlo solver.
//
// TODO: - Rename class.
//       - Const correctness.

#ifndef DCA_PHYS_MODELS_GENERAL_INTERACTION_HPP
#define DCA_PHYS_MODELS_GENERAL_INTERACTION_HPP

#include <array>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/convert.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename parameters_type>
class general_interaction {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;
  using r_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_dimension, domains::CLUSTER,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;

  template <class vertex_pair_type, class rng_type, class H_interaction_type>
  static void set_vertex(vertex_pair_type& vertex, parameters_type& parameters, rng_type& rng,
                         H_interaction_type H_interaction);

private:
  template <class H_interaction_type>
  static std::vector<int> make_correlated_orbitals(parameters_type& parameters,
                                                   H_interaction_type& H_interaction);

  static func::dmn_variadic<b, s, b, s, r_DCA, r_DCA>& get_domain() {
    static func::dmn_variadic<b, s, b, s, r_DCA, r_DCA> b_s_b_s_r_r_dmn;
    return b_s_b_s_r_r_dmn;
  }
};

template <typename parameters_type>
template <class vertex_pair_type, class rng_type, class H_interaction_type>
void general_interaction<parameters_type>::set_vertex(vertex_pair_type& vertex,
                                                      parameters_type& parameters, rng_type& rng,
                                                      H_interaction_type H_interaction) {
  static std::vector<int> correlated_orbitals =
      general_interaction<parameters_type>::make_correlated_orbitals(parameters, H_interaction);

  // Get a random pair of correlated orbitals
  const int pos = rng() * correlated_orbitals.size();
  const int lin_ind = correlated_orbitals[pos];

  std::array<int, 6> sub_ind;  // [0]=b1, [1]=s1, [2]=b2, [3]=s2, [4]=r1, [5]=r2
  get_domain().linind_2_subind(lin_ind, sub_ind.data());

  // Set the vertex properties
  vertex.get_bands().first = parameters.get_interacting_bands()[sub_ind[0]];
  vertex.get_bands().second = parameters.get_interacting_bands()[sub_ind[2]];

  vertex.get_e_spins().first = domains::electron_spin_domain::get_elements()[sub_ind[1]];
  vertex.get_e_spins().second = domains::electron_spin_domain::get_elements()[sub_ind[3]];

  vertex.get_spin_orbitals().first = domains::convert<int, nu>::spin_orbital(
      vertex.get_bands().first,
      vertex.get_e_spins().first);  // nu = func::dmn_variadic<b,s>
  vertex.get_spin_orbitals().second = domains::convert<int, nu>::spin_orbital(
      vertex.get_bands().second, vertex.get_e_spins().second);

  vertex.get_r_sites().first = sub_ind[4];
  vertex.get_r_sites().second = sub_ind[5];
}

template <typename parameters_type>
template <class H_interaction_type>
std::vector<int> general_interaction<parameters_type>::make_correlated_orbitals(
    parameters_type& /*parameters*/, H_interaction_type& H_interaction) {
  std::vector<int> correlated_orbitals;

  for (int r_j = 0; r_j < r_DCA::dmn_size(); ++r_j) {
    for (int r_i = 0; r_i < r_DCA::dmn_size(); ++r_i) {
      int delta_r = r_DCA::parameter_type::subtract(r_j, r_i);  // delta_r = r_i - r_j

      for (int s_j = 0; s_j < s::dmn_size(); ++s_j) {
        for (int b_j = 0; b_j < b::dmn_size(); ++b_j) {
          for (int s_i = 0; s_i < s::dmn_size(); ++s_i) {
            for (int b_i = 0; b_i < b::dmn_size(); ++b_i) {
              if (std::abs(H_interaction(b_i, s_i, b_j, s_j, delta_r)) > 1.e-3) {
                int linear_index = get_domain()(b_i, s_i, b_j, s_j, r_i, r_j);
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
