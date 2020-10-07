// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class kills old used HS-spins and resizes the N-matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/shrink_tools_algorithms/shrink_tools_algorithms.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, typename Real>
class SHRINK_TOOLS {
  typedef vertex_singleton vertex_singleton_type;
  using HostVector = linalg::util::HostVector<int>;

public:
  SHRINK_TOOLS(int id);

  template <class configuration_type, class vertex_vertex_matrix_type>
  static void shrink_Gamma(configuration_type& full_configuration,
                           vertex_vertex_matrix_type& Gamma_up, vertex_vertex_matrix_type& Gamma_dn);

  template <class configuration_type>
  void reorganize_configuration_test(configuration_type& full_configuration,
                                     dca::linalg::Matrix<Real, device_t>& N_up,
                                     dca::linalg::Matrix<Real, device_t>& N_dn,
                                     dca::linalg::Matrix<Real, device_t>& G0_up,
                                     dca::linalg::Matrix<Real, device_t>& G0_dn);

  int deviceFingerprint() const {
    return shrink_tools_algorithm_obj_.deviceFingerprint();
  }

private:
  template <class configuration_type, class vertex_vertex_matrix_type>
  static void shrink_Gamma_matrix(configuration_type& full_configuration,
                                  vertex_vertex_matrix_type& Gamma, e_spin_states_type e_spin);

  template <class configuration_type>
  static void shrink_Gamma_matrix(configuration_type& full_configuration,
                                  dca::linalg::Matrix<Real, device_t>& Gamma,
                                  e_spin_states_type e_spin);

  template <class configuration_type>
  static void swap_and_remove_vertices(configuration_type& full_configuration,
                                       linalg::util::HostVector<int>& source_index,
                                       linalg::util::HostVector<int>& target_index,
                                       e_spin_states_type e_spin);

  template <class configuration_type>
  static void erase_non_creatable_and_non_annihilatable_spins(
      configuration_type& full_configuration, dca::linalg::Matrix<Real, device_t>& N_up,
      dca::linalg::Matrix<Real, device_t>& N_dn, dca::linalg::Matrix<Real, device_t>& G0_up,
      dca::linalg::Matrix<Real, device_t>& G0_dn);

private:
  void test_swap_vectors(const linalg::util::HostVector<int>& source_index,
                         const linalg::util::HostVector<int>& target_index, int size);

private:
  int thread_id;
  int stream_id;

  HostVector source_index_up_;
  HostVector source_index_dn_;
  HostVector target_index_up_;
  HostVector target_index_dn_;

  SHRINK_TOOLS_ALGORITHMS<device_t, Real> shrink_tools_algorithm_obj_;
};

template <dca::linalg::DeviceType device_t, typename Real>
SHRINK_TOOLS<device_t, Real>::SHRINK_TOOLS(int id)
    : thread_id(id),
      stream_id(0),

      shrink_tools_algorithm_obj_(thread_id) {}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type, class vertex_vertex_matrix_type>
void SHRINK_TOOLS<device_t, Real>::shrink_Gamma(configuration_type& full_configuration,
                                                vertex_vertex_matrix_type& Gamma_up,
                                                vertex_vertex_matrix_type& Gamma_dn) {
  std::vector<int>& changed_spin_indices = full_configuration.get_changed_spin_indices();
  std::vector<HS_spin_states_type>& changed_spin_values =
      full_configuration.get_changed_spin_values();

  for (int i = 0; i < (int)changed_spin_indices.size();) {
    if (full_configuration[changed_spin_indices[i]].is_Bennett()) {
      changed_spin_indices.erase(changed_spin_indices.begin() + i);
      changed_spin_values.erase(changed_spin_values.begin() + i);
    }
    else
      i++;
  }

  shrink_Gamma_matrix(full_configuration, Gamma_up, e_UP);
  shrink_Gamma_matrix(full_configuration, Gamma_dn, e_DN);

  assert(2 * full_configuration.get_changed_spin_indices().size() ==
         full_configuration.get_changed_spin_indices_e_spin(e_UP).size() +
             full_configuration.get_changed_spin_indices_e_spin(e_DN).size());
}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type, class vertex_vertex_matrix_type>
void SHRINK_TOOLS<device_t, Real>::shrink_Gamma_matrix(configuration_type& full_configuration,
                                                       vertex_vertex_matrix_type& Gamma,
                                                       e_spin_states_type e_spin) {
  std::vector<int> configuration_spin_indices_to_be_erased(0);

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  std::vector<int>& changed_spin_indices_e_spin =
      full_configuration.get_changed_spin_indices_e_spin(e_spin);
  std::vector<HS_spin_states_type>& changed_spin_values_e_spin =
      full_configuration.get_changed_spin_values_e_spin(e_spin);

  for (int i = 0; i < Gamma.size();) {
    int configuration_index_e_spin = changed_spin_indices_e_spin[i];
    int configuration_index =
        configuration_e_spin[configuration_index_e_spin].get_configuration_index();

    if (full_configuration[configuration_index].is_Bennett()) {
      changed_spin_indices_e_spin.erase(changed_spin_indices_e_spin.begin() + i);
      changed_spin_values_e_spin.erase(changed_spin_values_e_spin.begin() + i);

      dca::linalg::matrixop::removeRowAndCol(Gamma, i);
    }
    else
      i++;
  }
}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type>
void SHRINK_TOOLS<device_t, Real>::shrink_Gamma_matrix(configuration_type& full_configuration,
                                                       dca::linalg::Matrix<Real, device_t>& Gamma,
                                                       e_spin_states_type e_spin) {
  std::vector<int> configuration_spin_indices_to_be_erased(0);

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  auto& changed_spin_indices_e_spin = full_configuration.get_changed_spin_indices_e_spin(e_spin);
  std::vector<HS_spin_states_type>& changed_spin_values_e_spin =
      full_configuration.get_changed_spin_values_e_spin(e_spin);

  assert(Gamma.is_square());

  for (int i = 0; i < Gamma.size().first;) {
    int configuration_index_e_spin = changed_spin_indices_e_spin[i];
    int configuration_index =
        configuration_e_spin[configuration_index_e_spin].get_configuration_index();

    if (full_configuration[configuration_index].is_Bennett()) {
      changed_spin_indices_e_spin.erase(changed_spin_indices_e_spin.begin() + i);
      changed_spin_values_e_spin.erase(changed_spin_values_e_spin.begin() + i);

      dca::linalg::matrixop::removeRowAndCol(Gamma, i);
    }
    else
      i++;
  }
}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type>
void SHRINK_TOOLS<device_t, Real>::reorganize_configuration_test(
    configuration_type& full_configuration, dca::linalg::Matrix<Real, device_t>& N_up,
    dca::linalg::Matrix<Real, device_t>& N_dn, dca::linalg::Matrix<Real, device_t>& G0_up,
    dca::linalg::Matrix<Real, device_t>& G0_dn) {
  source_index_up_.resize(0);
  target_index_up_.resize(0);
  source_index_dn_.resize(0);
  target_index_dn_.resize(0);

  swap_and_remove_vertices(full_configuration, source_index_up_, target_index_up_, e_UP);
  swap_and_remove_vertices(full_configuration, source_index_dn_, target_index_dn_, e_DN);

#ifndef NDEBUG
  test_swap_vectors(source_index_up_, target_index_up_, N_up.size().first);
  test_swap_vectors(source_index_dn_, target_index_dn_, N_dn.size().first);
#endif  // NDEBUG

  shrink_tools_algorithm_obj_.execute(source_index_up_, target_index_up_, N_up, G0_up,
                                      source_index_dn_, target_index_dn_, N_dn, G0_dn);

  erase_non_creatable_and_non_annihilatable_spins(full_configuration, N_up, N_dn, G0_up, G0_dn);
  assert(full_configuration.assert_consistency());
}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type>
void SHRINK_TOOLS<device_t, Real>::swap_and_remove_vertices(configuration_type& full_configuration,
                                                            HostVector& source_index,
                                                            HostVector& target_index,
                                                            e_spin_states_type e_spin) {
  const auto death_condition = [](const vertex_singleton& v) { return v.get_HS_spin() == HS_ZERO; };

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  if (configuration_size == 0) {
    return;
  }

  int dead_spin = 0;
  int living_spin = configuration_size - 1;

  auto update_config_idx = [&](int index, int new_index) {
    const int configuration_index = configuration_e_spin[index].get_configuration_index();
    auto& e_spin_indices = full_configuration[configuration_index].get_configuration_e_spin_indices();
    if (configuration_e_spin[index].get_HS_field() == HS_FIELD_DN) {
      e_spin_indices.first = new_index;
    }
    else {
      e_spin_indices.second = new_index;
    }
  };

  while (true) {
    while (dead_spin < configuration_size && !death_condition(configuration_e_spin[dead_spin])) {
      dead_spin++;
    }

    while (living_spin >= 0 && death_condition(configuration_e_spin[living_spin])) {
      update_config_idx(living_spin, -1);
      configuration_e_spin.pop_back();
      living_spin--;
    }

    if (dead_spin >= living_spin) {
      break;
    }

    assert(configuration_e_spin[dead_spin].get_HS_spin() == HS_ZERO);
    assert(configuration_e_spin[living_spin].get_HS_spin() != HS_ZERO);

    // Update index of vertex singletons in the vertex pair configuration.
    update_config_idx(living_spin, dead_spin);
    update_config_idx(dead_spin, -1);

    // Move living spin to the left and delete dead spin.
    configuration_e_spin[dead_spin] = configuration_e_spin[living_spin];
    configuration_e_spin.pop_back();

    // Store moves to apply to N and G0 matrices columns and rows.
    target_index.push_back(dead_spin);
    source_index.push_back(living_spin);

    dead_spin++;
    living_spin--;
  }
}

template <dca::linalg::DeviceType device_t, typename Real>
template <class configuration_type>
void SHRINK_TOOLS<device_t, Real>::erase_non_creatable_and_non_annihilatable_spins(
    configuration_type& full_configuration, dca::linalg::Matrix<Real, device_t>& N_up,
    dca::linalg::Matrix<Real, device_t>& N_dn, dca::linalg::Matrix<Real, device_t>& G0_up,
    dca::linalg::Matrix<Real, device_t>& G0_dn) {
  const int new_size_up = full_configuration.get(e_UP).size();
  N_up.resize(new_size_up);
  G0_up.resize(new_size_up);

  const int new_size_dn = full_configuration.get(e_DN).size();
  N_dn.resize(new_size_dn);
  G0_dn.resize(new_size_dn);

  for (unsigned i = 0; i < full_configuration.size();) {
    if (full_configuration[i].get_HS_spin() == HS_ZERO) {
      full_configuration.erase(i);
    }
    else {
      ++i;
    }
  }
}

template <dca::linalg::DeviceType device_t, typename Real>
void SHRINK_TOOLS<device_t, Real>::test_swap_vectors(const HostVector& source_index,
                                                     const HostVector& target_index, int size) {
  if (source_index.size() != target_index.size())
    throw std::logic_error("source_index.size() != target_index.size()");

  for (size_t i = 0; i < source_index.size(); ++i)
    if (source_index[i] < 0 or source_index[i] >= size)
      throw std::logic_error("source_index[i]<0 or source_index[i]>=size");

  for (size_t i = 0; i < target_index.size(); ++i)
    if (target_index[i] < 0 or target_index[i] >= size)
      throw std::logic_error("target_index[i]<0 or target_index[i]>=size");

  for (size_t i = 0; i < source_index.size(); ++i)
    for (size_t j = i + 1; j < source_index.size(); ++j)
      if (source_index[i] == source_index[j])
        throw std::logic_error("source_index[i] == source_index[j]");

  for (size_t i = 0; i < target_index.size(); ++i)
    for (size_t j = i + 1; j < target_index.size(); ++j)
      if (target_index[i] == target_index[j])
        throw std::logic_error("target_index[i] == target_index[j]");

  for (size_t i = 0; i < source_index.size(); ++i)
    for (size_t j = 0; j < target_index.size(); ++j)
      if (source_index[i] == target_index[j]) {
        for (size_t i = 0; i < source_index.size(); ++i)
          std::cout << i << "\t" << source_index[i] << "\t" << target_index[i] << std::endl;
        std::cout << std::endl;

        throw std::logic_error("source_index[i] == target_index[j]");
      }
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP
