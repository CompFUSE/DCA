// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the Hubbard Stratonovitch configuration.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CT_AUX_HS_CONFIGURATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CT_AUX_HS_CONFIGURATION_HPP

#include <cassert>
#include <cstdint>  // uint64_t
#include <cstdlib>  // std::size_t
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

#include "dca/io/buffer.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_field_sign_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_spin_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_pair.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class parameters_type>
class CT_AUX_HS_configuration {
public:
  using rng_type = typename parameters_type::random_number_generator;

  typedef HS_spin_states_type spin_state_type;

  typedef vertex_singleton vertex_singleton_type;
  typedef vertex_pair<parameters_type> vertex_pair_type;

public:
  CT_AUX_HS_configuration(parameters_type& parameters_ref, rng_type& rng_ref);

  int size();
  vertex_pair_type& operator[](int index);

  std::vector<vertex_singleton_type>& get(e_spin_states_type e_spin_type);
  const std::vector<vertex_singleton_type>& get(e_spin_states_type e_spin_type) const;

  void reset();

  // Creates an initial configuration with "initial-configuration-size" (input parameter) random
  // interacting vertices.
  void initialize();
  void prepare_configuration();
  void update_configuration_e_spin(vertex_pair_type& vertex_pair);

  void remove_HS_spin(int index);
  void add_delayed_HS_spin(int index, spin_state_type spin_value);
  void add_delayed_HS_spin_to_configuration_e_spin(int configuration_index,
                                                   HS_spin_states_type spin_value);

  void commit_accepted_spins();

  int get_first_non_interacting_spin_index(e_spin_states_type e_spin_type);
  int get_first_shuffled_spin_index(e_spin_states_type e_spin_type);

  std::vector<int>& get_changed_spin_indices();
  std::vector<HS_spin_states_type>& get_changed_spin_values();

  auto& get_changed_spin_indices_e_spin(e_spin_states_type e_spin_type);

  std::vector<HS_spin_states_type>& get_changed_spin_values_e_spin(e_spin_states_type e_spin_type);

  int get_number_of_interacting_HS_spins();

  int get_random_interacting_vertex();

  void insert_random_noninteracting_vertex(bool mark_annihilatable);

  // Remove element at index 'idx' in O(1). Vertices order is not preserved. The electron spin
  // configuration needs to be updated independently.
  void erase(unsigned idx);

  // debug tools
  void print() /*const*/;
  void print(e_spin_states_type e_spin_type) /*const*/;

  bool assert_block_form(
      e_spin_states_type e_spin_type) /*const*/;  // [non-shuffled-spin | shuffled-spins]
  bool assert_counters() /*const*/;
  bool assert_consistency() /*const*/;

  // Unmarks the vertex vertex_index as annihilatable.
  // Precondition: The vertex vertex_index is marked as annihilatable.
  // INTERNAL: This is a helper method to unmark all "virtual" interacting vertices, that have
  //           temporarily been marked as annihilatable by get_random_noninteracting_vertex.
  void unmarkAsAnnihilatable(const int vertex_index) {
    assert(configuration[vertex_index].is_annihilatable() == true);
    configuration[vertex_index].set_annihilatable(false);
    --current_Nb_of_annihilatable_spins_;
  }

  // Returns the position of the vertex with ID vertex_id or the size of the configuration if no
  // such vertex is found.
  std::size_t find(uint64_t vertex_id) const;

  bool operator==(const CT_AUX_HS_configuration<parameters_type>& rhs) const;

  template <class Pars>
  friend io::Buffer& operator<<(io::Buffer& buff, const CT_AUX_HS_configuration<Pars>& config);
  template <class Pars>
  friend io::Buffer& operator>>(io::Buffer& buff, CT_AUX_HS_configuration<Pars>& config);

private:
  parameters_type& parameters;
  rng_type& rng;

  std::vector<vertex_pair_type> configuration;

  std::vector<vertex_singleton_type> configuration_e_UP;  // = { configuration | e_spin == e_UP}
  std::vector<vertex_singleton_type> configuration_e_DN;  // = { configuration | e_spin == e_DN}

  unsigned current_Nb_of_annihilatable_spins_ = 0;

  std::vector<int> changed_spin_indices;
  std::vector<HS_spin_states_type> changed_spin_values;

  using HostVector = linalg::util::HostVector<int>;
  HostVector changed_spin_indices_e_UP;  // = { changed_spin_indices of configuration_e_UP}
  HostVector changed_spin_indices_e_DN;  // = { changed_spin_indices of configuration_e_DN}

  std::vector<HS_spin_states_type> changed_spin_values_e_UP;
  std::vector<HS_spin_states_type> changed_spin_values_e_DN;

  const int max_num_noninteracting_spins_;
  std::array<std::size_t, 2> size_at_last_step_;

  uint64_t next_vertex_id_;
};

template <class parameters_type>
CT_AUX_HS_configuration<parameters_type>::CT_AUX_HS_configuration(parameters_type& parameters_ref,
                                                                  rng_type& rng_ref)
    : parameters(parameters_ref),
      rng(rng_ref),

      // Rounding up ensures a value >= 1.
      max_num_noninteracting_spins_((parameters.get_max_submatrix_size() + 1) / 2),

      next_vertex_id_(0) {}

template <class parameters_type>
int CT_AUX_HS_configuration<parameters_type>::size() {
  return configuration.size();
}

template <class parameters_type>
typename CT_AUX_HS_configuration<parameters_type>::vertex_pair_type& CT_AUX_HS_configuration<
    parameters_type>::operator[](int index) {
  return configuration[index];
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::remove_HS_spin(int index) {
  // in this method, you assume that the configuration_e_UP/DN have been taken care of -->
  // SHRINK_TOOLS::shrink_configuration
  assert(!configuration[index].is_annihilatable() && !configuration[index].is_creatable());

  configuration.erase(configuration.begin() + index);

  for (size_t i = 0; i < configuration.size(); i++)
    if (configuration[i].get_configuration_index() > index)
      configuration[i].get_configuration_index() -= 1;

  for (size_t i = 0; i < configuration_e_UP.size(); i++)
    if (configuration_e_UP[i].get_configuration_index() > index)
      configuration_e_UP[i].get_configuration_index() -= 1;

  for (size_t i = 0; i < configuration_e_DN.size(); i++)
    if (configuration_e_DN[i].get_configuration_index() > index)
      configuration_e_DN[i].get_configuration_index() -= 1;
}

template <class parameters_type>
std::vector<vertex_singleton>& CT_AUX_HS_configuration<parameters_type>::get(e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return configuration_e_UP;
  else
    return configuration_e_DN;
}

template <class parameters_type>
const std::vector<vertex_singleton>& CT_AUX_HS_configuration<parameters_type>::get(
    e_spin_states_type e_spin) const {
  if (e_spin == e_UP)
    return configuration_e_UP;
  else
    return configuration_e_DN;
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::reset() {
  configuration.clear();

  configuration_e_UP.clear();
  configuration_e_DN.clear();

  current_Nb_of_annihilatable_spins_ = 0;

  changed_spin_indices.clear();
  changed_spin_values.clear();

  changed_spin_indices_e_UP.clear();
  changed_spin_values_e_UP.clear();

  changed_spin_indices_e_DN.clear();
  changed_spin_values_e_DN.clear();

  next_vertex_id_ = 0;
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::initialize() {
  reset();

  for (int i = 0; i < parameters.get_initial_configuration_size(); ++i) {
    vertex_pair_type vertex(parameters, rng, configuration.size(), next_vertex_id_++);
    vertex.set_random_interacting();

    ++current_Nb_of_annihilatable_spins_;
    configuration.push_back(vertex);
    update_configuration_e_spin(configuration.back());
  }
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::prepare_configuration() {
  assert(changed_spin_indices_e_UP.size() == 0);
  assert(changed_spin_indices_e_DN.size() == 0);
  assert(changed_spin_indices.size() == 0);

  size_at_last_step_ = {configuration_e_UP.size(), configuration_e_DN.size()};

  for (size_t i = 0; i < configuration.size(); i++) {
    configuration[i].set_shuffled(false);
    configuration[i].set_successfully_flipped(false);
    configuration[i].set_Bennett(false);

    assert(configuration[i].is_annihilatable());
  }

  assert(assert_consistency());
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::update_configuration_e_spin(
    vertex_pair_type& vertex_pair) {
  if (vertex_pair.get_e_spins().first == e_UP) {
    vertex_pair.get_configuration_e_spin_indices().first = configuration_e_UP.size();
    configuration_e_UP.push_back(vertex_pair.first());
  }
  else {
    vertex_pair.get_configuration_e_spin_indices().first = configuration_e_DN.size();
    configuration_e_DN.push_back(vertex_pair.first());
  }

  if (vertex_pair.get_e_spins().second == e_UP) {
    vertex_pair.get_configuration_e_spin_indices().second = configuration_e_UP.size();
    configuration_e_UP.push_back(vertex_pair.second());
  }
  else {
    vertex_pair.get_configuration_e_spin_indices().second = configuration_e_DN.size();
    configuration_e_DN.push_back(vertex_pair.second());
  }
}

template <class parameters_type>
int CT_AUX_HS_configuration<parameters_type>::get_first_non_interacting_spin_index(
    e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& spin_configuration = get(e_spin);

  // FIXME: What and when to return if configuration_size = 0?
  // Note: should be handled by the caller

  const int first_non_interacting_index = std::distance(
      spin_configuration.begin(),
      std::find_if(spin_configuration.begin(), spin_configuration.end(), [&](const auto& elem) {
        const int configuration_index = elem.get_configuration_index();
        return !configuration[configuration_index].is_annihilatable();
      }));

  return first_non_interacting_index;
}

template <class parameters_type>
int CT_AUX_HS_configuration<parameters_type>::get_first_shuffled_spin_index(e_spin_states_type e_spin) {
  // FIXME: What and when to return if configuration_size = 0?

  assert(assert_block_form(e_spin));

  //  const int first_shuffled_index = std::distance(
  //      spin_configuration.begin(),
  //      std::find_if(spin_configuration.begin(), spin_configuration.end(), [&](const auto& elem) {
  //        const int index = elem.get_configuration_index();
  //        return configuration[index].is_shuffled();
  //      }));

  return e_spin == e_UP ? size_at_last_step_[0] : size_at_last_step_[1];
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::commit_accepted_spins() {
  for (size_t i = 0; i < changed_spin_indices.size(); i++)
    configuration[changed_spin_indices[i]].get_HS_spin() = changed_spin_values[i];

  changed_spin_indices.clear();
  changed_spin_values.clear();

  for (size_t i = 0; i < changed_spin_indices_e_UP.size(); i++)
    configuration_e_UP[changed_spin_indices_e_UP[i]].get_HS_spin() = changed_spin_values_e_UP[i];

  changed_spin_indices_e_UP.clear();
  changed_spin_values_e_UP.clear();

  for (size_t i = 0; i < changed_spin_indices_e_DN.size(); i++)
    configuration_e_DN[changed_spin_indices_e_DN[i]].get_HS_spin() = changed_spin_values_e_DN[i];

  changed_spin_indices_e_DN.clear();
  changed_spin_values_e_DN.clear();
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::add_delayed_HS_spin(int configuration_index,
                                                                   spin_state_type spin_value) {
  add_delayed_HS_spin_to_configuration_e_spin(configuration_index, spin_value);

  if (configuration[configuration_index].is_successfully_flipped() &&
      spin_value == HS_ZERO)  // --> Bennett-case
  {
    for (size_t i = 0; i < changed_spin_indices.size(); i++)
      if (changed_spin_indices[i] == configuration_index)
        changed_spin_values[i] = spin_value;

    current_Nb_of_annihilatable_spins_ -= 1;

    configuration[configuration_index].set_annihilatable(false);
    configuration[configuration_index].set_Bennett(true);

    return;
  }

  changed_spin_indices.push_back(configuration_index);
  changed_spin_values.push_back(spin_value);

  if (spin_value == HS_ZERO) {
    // cout << "\t--> annihilate spin : " << configuration_index << std::endl;

    current_Nb_of_annihilatable_spins_ -= 1;

    configuration[configuration_index].set_annihilatable(false);
  }
  else {
    // cout << "\t--> create spin : " << configuration_index << std::endl;

    current_Nb_of_annihilatable_spins_ += 1;

    configuration[configuration_index].set_annihilatable(true);
  }

  configuration[configuration_index].set_successfully_flipped(true);
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::add_delayed_HS_spin_to_configuration_e_spin(
    int configuration_index, HS_spin_states_type spin_value) {
  auto& changed_spin_indices_e_spin_first =
      get_changed_spin_indices_e_spin(configuration[configuration_index].get_e_spins().first);
  std::vector<HS_spin_states_type>& changed_spin_values_e_spin_first =
      get_changed_spin_values_e_spin(configuration[configuration_index].get_e_spins().first);

  auto& changed_spin_indices_e_spin_second =
      get_changed_spin_indices_e_spin(configuration[configuration_index].get_e_spins().second);
  std::vector<HS_spin_states_type>& changed_spin_values_e_spin_second =
      get_changed_spin_values_e_spin(configuration[configuration_index].get_e_spins().second);

  int configuration_e_spin_index_first =
      configuration[configuration_index].get_configuration_e_spin_indices().first;
  int configuration_e_spin_index_second =
      configuration[configuration_index].get_configuration_e_spin_indices().second;

  if (configuration[configuration_index].is_successfully_flipped() &&
      spin_value == HS_ZERO)  // --> Bennett-case
  {
    for (size_t i = 0; i < changed_spin_indices_e_spin_first.size(); i++)
      if (changed_spin_indices_e_spin_first[i] == configuration_e_spin_index_first)
        changed_spin_values_e_spin_first[i] = spin_value;

    for (size_t i = 0; i < changed_spin_indices_e_spin_second.size(); i++)
      if (changed_spin_indices_e_spin_second[i] == configuration_e_spin_index_second)
        changed_spin_values_e_spin_second[i] = spin_value;

    return;  // -> do not push-back the Bennett spin !!
  }

  changed_spin_indices_e_spin_first.push_back(configuration_e_spin_index_first);
  changed_spin_values_e_spin_first.push_back(spin_value);

  changed_spin_indices_e_spin_second.push_back(configuration_e_spin_index_second);
  changed_spin_values_e_spin_second.push_back(spin_value);
}

template <class parameters_type>
std::vector<int>& CT_AUX_HS_configuration<parameters_type>::get_changed_spin_indices() {
  return changed_spin_indices;
}

template <class parameters_type>
std::vector<HS_spin_states_type>& CT_AUX_HS_configuration<parameters_type>::get_changed_spin_values() {
  return changed_spin_values;
}

template <class parameters_type>
auto& CT_AUX_HS_configuration<parameters_type>::get_changed_spin_indices_e_spin(
    e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return changed_spin_indices_e_UP;
  else
    return changed_spin_indices_e_DN;
}

template <class parameters_type>
std::vector<HS_spin_states_type>& CT_AUX_HS_configuration<parameters_type>::get_changed_spin_values_e_spin(
    e_spin_states_type e_spin) {
  if (e_spin == e_UP)
    return changed_spin_values_e_UP;
  else
    return changed_spin_values_e_DN;
}

template <class parameters_type>
int CT_AUX_HS_configuration<parameters_type>::get_random_interacting_vertex() {
  assert(current_Nb_of_annihilatable_spins_ > 0);

  int vertex_index;

  do {
    vertex_index = int(rng() * configuration.size());
  } while (!configuration[vertex_index].is_annihilatable());

  return vertex_index;
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::insert_random_noninteracting_vertex(
    bool mark_annihilatable) {
  vertex_pair_type vertex(parameters, rng, configuration.size(), next_vertex_id_++);
  vertex.set_random_noninteracting();

  configuration.push_back(vertex);
  update_configuration_e_spin(configuration.back());

  // However, this "virtual" interacting spin is eligible for removal.
  // INTERNAL: CtauxWalker::generateDelayedSpinsAbortAtBennett unmarks the "virtual" interacting
  //           spins as annihilatable when all delayed spins have been generated.
  if (mark_annihilatable) {
    configuration.back().set_annihilatable(true);
    ++current_Nb_of_annihilatable_spins_;
  }
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::erase(unsigned idx) {
  assert(configuration[idx].get_HS_spin() == HS_ZERO);
  assert(configuration[idx].get_configuration_e_spin_indices().first == -1);
  assert(configuration[idx].get_configuration_e_spin_indices().second == -1);

  if (configuration[idx].is_annihilatable())
    --current_Nb_of_annihilatable_spins_;

  configuration[idx] = configuration.back();
  configuration.pop_back();
  if (idx == size())
    return;

  const auto& e_spin_indices = configuration[idx].get_configuration_e_spin_indices();
  const auto& e_spins = configuration[idx].get_e_spins();

  const auto update_index = [&](const auto& e_spin, const int spin_index) {
    if (spin_index >= 0) {
      assert(spin_index < get(e_spin).size());
      get(e_spin)[spin_index].get_configuration_index() = idx;
    }
  };

  update_index(e_spins.first, e_spin_indices.first);
  update_index(e_spins.second, e_spin_indices.second);

  configuration[idx].get_configuration_index() = idx;
}

template <class parameters_type>
int CT_AUX_HS_configuration<parameters_type>::get_number_of_interacting_HS_spins() {
  assert(assert_counters());
  return current_Nb_of_annihilatable_spins_;
}

template <class parameters_type>
bool CT_AUX_HS_configuration<parameters_type>::assert_block_form(e_spin_states_type e_spin) /*const*/ {
  const std::vector<vertex_singleton_type>& configuration_e_spin = get(e_spin);
  int configuration_size = configuration_e_spin.size();

  int vertex_index = 0;
  while (vertex_index < configuration_size &&
         configuration_e_spin[vertex_index].get_HS_spin() != HS_ZERO)
    vertex_index++;

  while (vertex_index < configuration_size &&
         configuration_e_spin[vertex_index].get_HS_spin() == HS_ZERO)
    vertex_index++;

  if (vertex_index != configuration_size) {
    std::cout << vertex_index << " --> " << configuration_size << std::endl;
    print();
    print(e_spin);
  }
  assert(vertex_index == configuration_size);

  return true;
}

template <class parameters_type>
bool CT_AUX_HS_configuration<parameters_type>::assert_counters() /*const*/ {
  int tmp = 0;

  for (size_t i = 0; i < configuration.size(); i++) {
    if (configuration[i].is_annihilatable())
      ++tmp;
  }

  if (tmp != current_Nb_of_annihilatable_spins_)
    throw std::logic_error("tmp != current_Nb_of_annihilatable_spins");
  return true;
}

template <class parameters_type>
bool CT_AUX_HS_configuration<parameters_type>::assert_consistency() /*const*/ {
  assert(2 * configuration.size() == (configuration_e_UP.size() + configuration_e_DN.size()));

  assert_counters();

  // assert configuration-index
  for (int i = 0; i < (int)configuration.size(); i++)
    assert(configuration[i].get_configuration_index() == i);

  {  // assert get_configuration_e_spin_indices() && get_configuration_index()

    for (int i = 0; i < (int)configuration_e_UP.size(); i++) {
      int configuration_index = configuration_e_UP[i].get_configuration_index();
      HS_field_sign_type HS_field = configuration_e_UP[i].get_HS_field();

      if (HS_field == HS_FIELD_DN) {
        // assert(configuration[configuration_index].get_configuration_e_spin_indices().first  ==
        // i);
        if (configuration[configuration_index].get_configuration_e_spin_indices().first != i)
          throw std::logic_error(__FUNCTION__);
      }
      else {
        // assert(configuration[configuration_index].get_configuration_e_spin_indices().second ==
        // i);
        if (configuration[configuration_index].get_configuration_e_spin_indices().second != i)
          throw std::logic_error(__FUNCTION__);
      }
    }

    for (int i = 0; i < (int)configuration_e_DN.size(); i++) {
      int configuration_index = configuration_e_DN[i].get_configuration_index();
      HS_field_sign_type HS_field = configuration_e_DN[i].get_HS_field();

      if (HS_field == HS_FIELD_DN) {
        // assert(configuration[configuration_index].get_configuration_e_spin_indices().first  ==
        // i);
        if (configuration[configuration_index].get_configuration_e_spin_indices().first != i)
          throw std::logic_error(__FUNCTION__);
      }
      else {
        // assert(configuration[configuration_index].get_configuration_e_spin_indices().second ==
        // i);
        if (configuration[configuration_index].get_configuration_e_spin_indices().second != i)
          throw std::logic_error(__FUNCTION__);
      }
    }
  }

  {  // assert internal pointers
    for (int i = 0; i < (int)configuration_e_UP.size(); i++) {
      const vertex_singleton_type& partner = configuration_e_UP[i].get_partner(*this);
      // assert( partner.get_configuration_index() ==
      // configuration_e_UP[i].get_configuration_index());
      if (partner.get_configuration_index() != configuration_e_UP[i].get_configuration_index())
        throw std::logic_error(__FUNCTION__);
    }

    for (int i = 0; i < (int)configuration_e_DN.size(); i++) {
      const vertex_singleton_type& partner = configuration_e_DN[i].get_partner(*this);
      // assert( partner.get_configuration_index() ==
      // configuration_e_DN[i].get_configuration_index());
      if (partner.get_configuration_index() != configuration_e_DN[i].get_configuration_index())
        throw std::logic_error(__FUNCTION__);
    }
  }

  return true;
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::print() /*const*/ {
  std::stringstream ss;
  ss << std::scientific;
  ss.precision(6);
  ss.width(6);

  ss << "\n         ";
  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t"
       << "==============";

  std::cout << "current_Nb_of_annihilatable_spins \t" << current_Nb_of_annihilatable_spins_
            << std::endl;
  std::cout << std::endl;

  ss << "\n tau     ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t" << configuration[i].get_tau();

  ss << "\n HS_spin ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t\t" << configuration[i].get_HS_spin();

  ss << "\n sites   ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t < " << configuration[i].get_r_sites().first << " , "
       << configuration[i].get_r_sites().second << " >";

  ss << "\n s-o     ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t < " << configuration[i].get_spin_orbitals().first << " , "
       << configuration[i].get_spin_orbitals().second << " >";

  ss << "\n bands   ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t < " << configuration[i].get_bands().first << " , "
       << configuration[i].get_bands().second << " >";

  ss << "\n es      ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t < " << configuration[i].get_e_spins().first << " , "
       << configuration[i].get_e_spins().second << " >";

  ss << "\n indices ";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t < " << configuration[i].get_configuration_e_spin_indices().first << " , "
       << configuration[i].get_configuration_e_spin_indices().second << " >";

  ss << "\n         ";
  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t"
       << "--------------";

  ss << "\n annihilatable\n\t";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t\t" << configuration[i].is_annihilatable();

  ss << "\n";

  ss << "\n changed\n\t";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t\t" << configuration[i].is_successfully_flipped();

  ss << "\n";

  ss << "\n Bennett\n\t";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t\t" << configuration[i].is_Bennett();

  ss << "\n";

  ss << "\n shuffled\n\t";

  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t\t" << configuration[i].is_shuffled();

  ss << "\n         ";
  for (size_t i = 0; i < configuration.size(); i++)
    ss << "\t"
       << "==============";

  ss << "\n";

  ss << std::endl << std::endl;

  std::cout << ss.str();
}

template <class parameters_type>
void CT_AUX_HS_configuration<parameters_type>::print(e_spin_states_type e_spin) /*const*/ {
  std::stringstream ss;
  ss << std::scientific;
  ss.precision(6);
  ss.width(6);

  ss << "\t  e_spin_states_type :: " << e_spin;

  ss << "\n               ";
  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t"
       << "==============";

  ss << "\n tau           ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t" << get(e_spin)[i].get_tau();

  ss << "\n site          ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_r_site();

  ss << "\n s-o           ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_spin_orbital();

  ss << "\n b             ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_band();

  ss << "\n e_spin        ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_e_spin();

  ss << "\n HS_spin       ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_HS_spin();

  ss << "\n HS_field_sign ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_HS_field();

  ss << "\n index         ";

  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t\t" << get(e_spin)[i].get_configuration_index();

  ss << "\n         ";
  for (size_t i = 0; i < get(e_spin).size(); i++)
    ss << "\t"
       << "==============";

  ss << std::endl << std::endl;

  std::cout << ss.str();
}

template <class parameters_type>
std::size_t CT_AUX_HS_configuration<parameters_type>::find(const uint64_t vertex_id) const {
  auto it = std::find_if(
      configuration.begin(), configuration.end(),
      [vertex_id](const vertex_pair_type& vertex) -> bool { return vertex.get_id() == vertex_id; });

  return it - configuration.begin();
}

template <class parameters_type>
bool CT_AUX_HS_configuration<parameters_type>::operator==(
    const CT_AUX_HS_configuration<parameters_type>& rhs) const {
  using math::util::operator==;

  // assert(assert_counters());
  return configuration == rhs.configuration && configuration_e_UP == rhs.configuration_e_UP &&
         configuration_e_DN == rhs.configuration_e_DN;
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CT_AUX_HS_CONFIGURATION_HPP
