// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a vertex pair consisting of two vertex singlets.
//
// TODO: Make this class fulfill the rule of 3 (5).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_PAIR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_PAIR_HPP

#include <cassert>
#include <cstdint>  // uint64_t
#include <utility>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/ctaux/structs/cv.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/models/general_interaction.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class parameters_type>
class vertex_pair {
public:
  using rng_type = typename parameters_type::random_number_generator;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  typedef RClusterDmn r_dmn_t;
  typedef typename r_dmn_t::parameter_type r_cluster_type;

  typedef vertex_singleton vertex_singleton_type;
  typedef vertex_pair<parameters_type> this_type;

public:
  vertex_pair(parameters_type& parameters_ref, rng_type& rng_ref, int configuration_index_in,
              int configuration_e_DN_index_in, int configuration_e_UP_index_in, uint64_t id);

  this_type& operator=(const this_type& other_vertex_pair);

  uint64_t get_id() const {
    return id_;
  }

  vertex_singleton_type first();
  vertex_singleton_type second();

  void set_random_interacting();
  void set_random_noninteracting();

  // physics-information
  std::pair<int, int>& get_bands();
  std::pair<e_spin_states_type, e_spin_states_type>& get_e_spins();
  std::pair<int, int>& get_r_sites();
  std::pair<int, int>& get_spin_orbitals();

  HS_spin_states_type& get_HS_spin();
  HS_spin_states_type get_HS_spin() const;
  int& get_delta_r();
  double& get_tau();

  // algorithm-information

  int& get_configuration_index();
  std::pair<int, int>& get_configuration_e_spin_indices();

  bool& is_creatable();
  bool& is_annihilatable();
  bool& is_shuffled();
  bool& is_Bennett();
  bool& is_successfully_flipped();

  bool operator==(const vertex_pair<parameters_type>& rhs) const;

  template <class Pars>
  friend io::Buffer& operator<<(io::Buffer& buff, const vertex_pair<Pars>& config);
  template <class Pars>
  friend io::Buffer& operator>>(io::Buffer& buff, vertex_pair<Pars>& config);

private:
  parameters_type& parameters;
  //     concurrency_type&   concurrency;
  rng_type& rng;

  uint64_t id_;

  const std::vector<int>& interacting_bands;
  int BANDS;

  // physics-information
  std::pair<int, int> bands;
  std::pair<e_spin_states_type, e_spin_states_type> e_spins;
  std::pair<int, int> spin_orbitals;
  std::pair<HS_field_sign_type, HS_field_sign_type> HS_fields;
  std::pair<int, int> r_sites;

  HS_spin_states_type HS_spin;
  int delta_r;
  double tau;

  // algorithm-information
  int configuration_index;
  std::pair<int, int> configuration_e_spin_indices;

  bool creatable;
  bool annihilatable;
  bool successfully_flipped;
  bool Bennett;
  bool shuffled;
};

template <class parameters_type>
vertex_pair<parameters_type>::vertex_pair(parameters_type& parameters_ref, rng_type& rng_ref,
                                          int configuration_index_in,
                                          int /*configuration_e_DN_index_in*/,
                                          int /*configuration_e_UP_index_in*/, uint64_t id)
    : parameters(parameters_ref),
      rng(rng_ref),
      //     concurrency(parameters.get_concurrency()),

      id_(id),

      interacting_bands(parameters_ref.get_interacting_orbitals()),
      BANDS(interacting_bands.size()),

      bands(0, 0),
      e_spins(e_DN, e_UP),
      spin_orbitals(0, 1),
      HS_fields(HS_FIELD_DN, HS_FIELD_UP),
      r_sites(0, 0),

      HS_spin(HS_ZERO),
      delta_r(0),
      tau(0),

      configuration_index(configuration_index_in),

      configuration_e_spin_indices(-1, -1),

      creatable(true),
      annihilatable(false),
      successfully_flipped(false),
      Bennett(false),
      shuffled(false) {}

template <class parameters_type>
vertex_pair<parameters_type>& vertex_pair<parameters_type>::operator=(
    const vertex_pair<parameters_type>& other_vertex_pair)  // --> necessary for push_back
{
  id_ = other_vertex_pair.id_;

  bands = other_vertex_pair.bands;
  e_spins = other_vertex_pair.e_spins;
  spin_orbitals = other_vertex_pair.spin_orbitals;
  r_sites = other_vertex_pair.r_sites;

  HS_spin = other_vertex_pair.HS_spin;
  delta_r = other_vertex_pair.delta_r;
  tau = other_vertex_pair.tau;

  configuration_index = other_vertex_pair.configuration_index;
  configuration_e_spin_indices = other_vertex_pair.configuration_e_spin_indices;

  creatable = other_vertex_pair.creatable;
  annihilatable = other_vertex_pair.annihilatable;
  successfully_flipped = other_vertex_pair.successfully_flipped;
  Bennett = other_vertex_pair.Bennett;
  shuffled = other_vertex_pair.shuffled;

  return *this;
}

template <class parameters_type>
vertex_singleton vertex_pair<parameters_type>::first() {
  return vertex_singleton_type(bands.first, e_spins.first, spin_orbitals.first,
                               spin_orbitals.second, r_sites.first, delta_r, tau, HS_spin,
                               HS_FIELD_DN, configuration_index);
}

template <class parameters_type>
vertex_singleton vertex_pair<parameters_type>::second() {
  return vertex_singleton_type(bands.second, e_spins.second, spin_orbitals.second,
                               spin_orbitals.first, r_sites.second, delta_r, tau, HS_spin,
                               HS_FIELD_UP, configuration_index);
}

template <class parameters_type>
void vertex_pair<parameters_type>::set_random_interacting() {
  dca::phys::models::general_interaction<parameters_type>::set_vertex(
      *this, parameters, rng, CV<parameters_type>::get_H_interaction());

  double draw = rng();

  if (draw > 1 / 2.)
    HS_spin = HS_UP;
  else
    HS_spin = HS_DN;

  delta_r = r_cluster_type::subtract(r_sites.second, r_sites.first);  // delta_r = r_i - r_j

  tau = parameters.get_beta() * rng();

  creatable = false;
  annihilatable = true;

  successfully_flipped = false;
  Bennett = false;
  shuffled = true;

  // For on-site interaction the spin-orbitals must be different.
  assert(delta_r != 0 || spin_orbitals.first != spin_orbitals.second);
}

template <class parameters_type>
void vertex_pair<parameters_type>::set_random_noninteracting() {
  dca::phys::models::general_interaction<parameters_type>::set_vertex(
      *this, parameters, rng, CV<parameters_type>::get_H_interaction());

  HS_spin = HS_ZERO;

  delta_r = r_cluster_type::subtract(r_sites.second, r_sites.first);  // delta_r = r_i - r_j

  tau = parameters.get_beta() * rng();

  creatable = true;
  annihilatable = false;

  successfully_flipped = false;
  Bennett = false;
  shuffled = true;

  assert(delta_r != 0 || spin_orbitals.first != spin_orbitals.second);
}

template <class parameters_type>
std::pair<int, int>& vertex_pair<parameters_type>::get_bands() {
  return bands;
}

template <class parameters_type>
std::pair<e_spin_states_type, e_spin_states_type>& vertex_pair<parameters_type>::get_e_spins() {
  return e_spins;
}

template <class parameters_type>
std::pair<int, int>& vertex_pair<parameters_type>::get_spin_orbitals() {
  return spin_orbitals;
}

template <class parameters_type>
std::pair<int, int>& vertex_pair<parameters_type>::get_r_sites() {
  return r_sites;
}

template <class parameters_type>
HS_spin_states_type& vertex_pair<parameters_type>::get_HS_spin() {
  return HS_spin;
}

template <class parameters_type>
HS_spin_states_type vertex_pair<parameters_type>::get_HS_spin() const {
  return HS_spin;
}

template <class parameters_type>
int& vertex_pair<parameters_type>::get_delta_r() {
  return delta_r;
}

template <class parameters_type>
double& vertex_pair<parameters_type>::get_tau() {
  return tau;
}

template <class parameters_type>
int& vertex_pair<parameters_type>::get_configuration_index() {
  return configuration_index;
}

template <class parameters_type>
std::pair<int, int>& vertex_pair<parameters_type>::get_configuration_e_spin_indices() {
  return configuration_e_spin_indices;
}

template <class parameters_type>
bool& vertex_pair<parameters_type>::is_creatable() {
  return creatable;
}

template <class parameters_type>
bool& vertex_pair<parameters_type>::is_annihilatable() {
  return annihilatable;
}

template <class parameters_type>
bool& vertex_pair<parameters_type>::is_successfully_flipped() {
  return successfully_flipped;
}

template <class parameters_type>
bool& vertex_pair<parameters_type>::is_Bennett() {
  return Bennett;
}

template <class parameters_type>
bool& vertex_pair<parameters_type>::is_shuffled() {
  return shuffled;
}

template <class parameters_type>
bool vertex_pair<parameters_type>::operator==(
    const dca::phys::solver::ctaux::vertex_pair<parameters_type>& rhs) const {
  using math::util::operator==;

  return id_ == rhs.id_ && interacting_bands == rhs.interacting_bands && BANDS == rhs.BANDS &&
         bands == rhs.bands && e_spins == rhs.e_spins && spin_orbitals == rhs.spin_orbitals &&
         HS_fields == rhs.HS_fields && r_sites == rhs.r_sites && HS_spin == rhs.HS_spin &&
         delta_r == rhs.delta_r && tau == rhs.tau && configuration_index == rhs.configuration_index &&
         configuration_e_spin_indices == rhs.configuration_e_spin_indices &&
         creatable == rhs.creatable && annihilatable == rhs.annihilatable &&
         successfully_flipped == rhs.successfully_flipped && Bennett == rhs.Bennett &&
         shuffled == rhs.shuffled;
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_PAIR_HPP
