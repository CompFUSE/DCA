// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class kills old used HS-spins and resizes the N-matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/shrink_tools_algorithms/shrink_tools_algorithms.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t>
class SHRINK_TOOLS {
  typedef vertex_singleton vertex_singleton_type;

public:
  SHRINK_TOOLS(int id);

  template <class configuration_type, class vertex_vertex_matrix_type>
  static void shrink_Gamma(configuration_type& full_configuration,
                           vertex_vertex_matrix_type& Gamma_up, vertex_vertex_matrix_type& Gamma_dn);

  /*
    template<class configuration_type, class vertex_vertex_matrix_type>
    static void reorganize_configuration(configuration_type&        full_configuration,
    vertex_vertex_matrix_type& N_up,
    vertex_vertex_matrix_type& N_dn,
    vertex_vertex_matrix_type& G0_up,
    vertex_vertex_matrix_type& G0_dn);
  */

  template <class configuration_type>
  void reorganize_configuration_test(configuration_type& full_configuration,
                                     dca::linalg::Matrix<double, device_t>& N_up,
                                     dca::linalg::Matrix<double, device_t>& N_dn,
                                     dca::linalg::Matrix<double, device_t>& G0_up,
                                     dca::linalg::Matrix<double, device_t>& G0_dn);

  int deviceFingerprint() const {
    return SHRINK_TOOLS_ALGORITHMS_obj.deviceFingerprint();
  }

private:
  template <class configuration_type, class vertex_vertex_matrix_type>
  static void shrink_Gamma_matrix(configuration_type& full_configuration,
                                  vertex_vertex_matrix_type& Gamma, e_spin_states_type e_spin);

  template <class configuration_type>
  static void shrink_Gamma_matrix(configuration_type& full_configuration,
                                  dca::linalg::Matrix<double, device_t>& Gamma,
                                  e_spin_states_type e_spin);

  /*
    template<class configuration_type>
    void swap_interacting_vertices_to_left(configuration_type&        full_configuration,
    dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G0,
    e_spin_states_type         e_spin);
  */

  template <class configuration_type>
  static void swap_interacting_vertices_to_left(configuration_type& full_configuration,
                                                std::vector<int>& source_index,
                                                std::vector<int>& target_index,
                                                e_spin_states_type e_spin);

  /*
    template<class configuration_type>
    static void swap_non_changed_vertices_to_left(configuration_type&        full_configuration,
    dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G0,
    e_spin_states_type         e_spin);
  */

  template <class configuration_type>
  static void swap_non_changed_vertices_to_left(configuration_type& full_configuration,
                                                std::vector<int>& source_index,
                                                std::vector<int>& target_index,
                                                e_spin_states_type e_spin);

  template <class configuration_type>
  static void erase_non_creatable_and_non_annihilatable_spins(
      configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N_up,
      dca::linalg::Matrix<double, device_t>& N_dn, dca::linalg::Matrix<double, device_t>& G0_up,
      dca::linalg::Matrix<double, device_t>& G0_dn);

  template <class configuration_type>
  static void erase_non_creatable_and_non_annihilatable_spins(
      configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N,
      dca::linalg::Matrix<double, device_t>& G0, e_spin_states_type e_spin);

private:
  bool test_swap_vectors(std::vector<int>& source_index, std::vector<int>& target_index, int size);

private:
  int thread_id;
  int stream_id;

  std::vector<int> source_index_up;
  std::vector<int> source_index_dn;

  std::vector<int> target_index_up;
  std::vector<int> target_index_dn;

  SHRINK_TOOLS_ALGORITHMS<device_t> SHRINK_TOOLS_ALGORITHMS_obj;
};

template <dca::linalg::DeviceType device_t>
SHRINK_TOOLS<device_t>::SHRINK_TOOLS(int id)
    : thread_id(id),
      stream_id(0),

      source_index_up(0),
      source_index_dn(0),

      target_index_up(0),
      target_index_dn(0),

      SHRINK_TOOLS_ALGORITHMS_obj(thread_id) {}

template <dca::linalg::DeviceType device_t>
template <class configuration_type, class vertex_vertex_matrix_type>
void SHRINK_TOOLS<device_t>::shrink_Gamma(configuration_type& full_configuration,
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

template <dca::linalg::DeviceType device_t>
template <class configuration_type, class vertex_vertex_matrix_type>
void SHRINK_TOOLS<device_t>::shrink_Gamma_matrix(configuration_type& full_configuration,
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

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::shrink_Gamma_matrix(configuration_type& full_configuration,
                                                 dca::linalg::Matrix<double, device_t>& Gamma,
                                                 e_spin_states_type e_spin) {
  std::vector<int> configuration_spin_indices_to_be_erased(0);

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  std::vector<int>& changed_spin_indices_e_spin =
      full_configuration.get_changed_spin_indices_e_spin(e_spin);
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

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::reorganize_configuration_test(
    configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N_up,
    dca::linalg::Matrix<double, device_t>& N_dn, dca::linalg::Matrix<double, device_t>& G0_up,
    dca::linalg::Matrix<double, device_t>& G0_dn) {
  {
    source_index_up.resize(0);
    target_index_up.resize(0);

    source_index_dn.resize(0);
    target_index_dn.resize(0);

    swap_interacting_vertices_to_left(full_configuration, source_index_up, target_index_up, e_UP);
    swap_interacting_vertices_to_left(full_configuration, source_index_dn, target_index_dn, e_DN);

    // cout << "swap_interacting_vertices_to_left" << endl;
    test_swap_vectors(source_index_up, target_index_up, N_up.size().first);
    test_swap_vectors(source_index_dn, target_index_dn, N_dn.size().first);

    if (false) {
      // SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_up, target_index_up, N_up, G0_up);
      // SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_dn, target_index_dn, N_dn, G0_dn);
    }
    else {
      SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_up, target_index_up, N_up, G0_up,
                                          source_index_dn, target_index_dn, N_dn, G0_dn);
    }
  }

  {
    source_index_up.resize(0);
    target_index_up.resize(0);

    source_index_dn.resize(0);
    target_index_dn.resize(0);

    swap_non_changed_vertices_to_left(full_configuration, source_index_up, target_index_up, e_UP);
    swap_non_changed_vertices_to_left(full_configuration, source_index_dn, target_index_dn, e_DN);

    // cout << "swap_non_changed_vertices_to_left" << endl;
    test_swap_vectors(source_index_up, target_index_up, N_up.size().first);
    test_swap_vectors(source_index_dn, target_index_dn, N_dn.size().first);

    if (false) {
      // SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_up, target_index_up, N_up, G0_up);
      // SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_dn, target_index_dn, N_dn, G0_dn);
    }
    else {
      SHRINK_TOOLS_ALGORITHMS_obj.execute(source_index_up, target_index_up, N_up, G0_up,
                                          source_index_dn, target_index_dn, N_dn, G0_dn);
    }
  }

  erase_non_creatable_and_non_annihilatable_spins(full_configuration, N_up, N_dn, G0_up, G0_dn);

  assert(full_configuration.assert_consistency());

  /*
    for(size_t l=0; l<source_index_up.size(); ++l){
    dca::linalg::matrixop::swapRowAndCol(N_up , source_index_up[l], target_index_up[l]);
    dca::linalg::matrixop::swapRowAndCol(G0_up, source_index_up[l], target_index_up[l]);
    }

    for(size_t l=0; l<source_index_dn.size(); ++l){
    dca::linalg::matrixop::swapRowAndCol(N_dn , source_index_dn[l], target_index_dn[l]);
    dca::linalg::matrixop::swapRowAndCol(G0_dn, source_index_dn[l], target_index_dn[l]);
    }
  */
}

/*
  template<dca::linalg::DeviceType device_t>
  template<class configuration_type>
  void SHRINK_TOOLS<device_t>::swap_interacting_vertices_to_left(configuration_type&
  full_configuration,
  dca::linalg::Matrix<double, device_t>& N,
  dca::linalg::Matrix<double, device_t>& G0,
  e_spin_states_type                 e_spin)
  {
  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int                                 configuration_size   = configuration_e_spin.size();

  int dead_spin   = 0;
  int configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();

  int living_spin = configuration_size-1;
  int configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();

  while(true)
  {
  while(dead_spin < configuration_size-1
  &&  full_configuration[configuration_index_dead_spin].is_annihilatable() )
  {
  dead_spin++;
  configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
  }

  while(living_spin > 0
  && !full_configuration[configuration_index_living_spin].is_annihilatable() )
  {
  living_spin--;
  configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
  }

  if(dead_spin > living_spin)
  break;
  else
  {
  HS_field_sign_type HS_field_dead_spin   = configuration_e_spin[dead_spin  ].get_HS_field();
  HS_field_sign_type HS_field_living_spin = configuration_e_spin[living_spin].get_HS_field();

  std::pair<int,int>& pair_dead_spin   = full_configuration[configuration_index_dead_spin]
  .get_configuration_e_spin_indices();
  std::pair<int,int>& pair_living_spin =
  full_configuration[configuration_index_living_spin].get_configuration_e_spin_indices();

  dca::linalg::matrixop::swapRowAndCol(N , dead_spin, living_spin, thread_id, stream_id);
  dca::linalg::matrixop::swapRowAndCol(G0, dead_spin, living_spin, thread_id, stream_id);

  swap(configuration_e_spin[dead_spin], configuration_e_spin[living_spin]);

  if(HS_field_dead_spin == HS_FIELD_DN)
  pair_dead_spin.first = living_spin;
  else
  pair_dead_spin.second = living_spin;

  if(HS_field_living_spin == HS_FIELD_DN)
  pair_living_spin.first = dead_spin;
  else
  pair_living_spin.second = dead_spin;

  }

  dead_spin++;
  living_spin--;

  if(dead_spin == int(configuration_e_spin.size()) || living_spin == -1)
  break;
  else
  {
  configuration_index_dead_spin   = configuration_e_spin[dead_spin]  .get_configuration_index();
  configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
  }
  }
  }
*/

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::swap_interacting_vertices_to_left(configuration_type& full_configuration,
                                                               std::vector<int>& source_index,
                                                               std::vector<int>& target_index,
                                                               e_spin_states_type e_spin) {
  // cout << __FUNCTION__ << endl;

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  if (configuration_size == 0) {
    return;
  }

  int dead_spin = 0;
  int configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();

  int living_spin = configuration_size - 1;
  int configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();

  while (true) {
    while (dead_spin < configuration_size - 1 &&
           full_configuration[configuration_index_dead_spin].is_annihilatable()) {
      dead_spin++;
      configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
    }

    while (living_spin > 0 &&
           !full_configuration[configuration_index_living_spin].is_annihilatable()) {
      living_spin--;
      configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
    }

    // if(dead_spin > living_spin)
    if (dead_spin >= living_spin)
      break;
    else {
      HS_field_sign_type HS_field_dead_spin = configuration_e_spin[dead_spin].get_HS_field();
      HS_field_sign_type HS_field_living_spin = configuration_e_spin[living_spin].get_HS_field();

      std::pair<int, int>& pair_dead_spin =
          full_configuration[configuration_index_dead_spin].get_configuration_e_spin_indices();
      std::pair<int, int>& pair_living_spin =
          full_configuration[configuration_index_living_spin].get_configuration_e_spin_indices();

      source_index.push_back(dead_spin);
      target_index.push_back(living_spin);

      // dca::linalg::matrixop::swapRowAndCol(N , dead_spin, living_spin);
      // dca::linalg::matrixop::swapRowAndCol(G0, dead_spin, living_spin);

      std::swap(configuration_e_spin[dead_spin], configuration_e_spin[living_spin]);

      if (HS_field_dead_spin == HS_FIELD_DN)
        pair_dead_spin.first = living_spin;
      else
        pair_dead_spin.second = living_spin;

      if (HS_field_living_spin == HS_FIELD_DN)
        pair_living_spin.first = dead_spin;
      else
        pair_living_spin.second = dead_spin;
    }

    dead_spin++;
    living_spin--;

    if (dead_spin == int(configuration_e_spin.size()) || living_spin == -1)
      break;
    else {
      configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
      configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
    }
  }
}

/*
  template<dca::linalg::DeviceType device_t>
  template<class configuration_type>
  void SHRINK_TOOLS<device_t>::swap_non_changed_vertices_to_left(configuration_type&
  full_configuration,
  dca::linalg::Matrix<double, device_t>& N,
  dca::linalg::Matrix<double, device_t>& G0,
  e_spin_states_type         e_spin)
  {
  //cout << __FUNCTION__ << endl;

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int                            configuration_size   = configuration_e_spin.size();

  int dead_spin   = 0;
  int configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();

  int living_spin = configuration_size-1;
  int configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();

  while(true)
  {
  while( dead_spin < configuration_size-1
  &&!(!full_configuration[configuration_index_dead_spin].is_annihilatable()
  && !full_configuration[configuration_index_dead_spin].is_creatable()) )
  {
  dead_spin++;
  configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
  }

  while(living_spin > 0
  && !full_configuration[configuration_index_living_spin].is_creatable() )
  {
  living_spin--;
  configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
  }

  if(dead_spin > living_spin)
  break;
  else
  {
  HS_field_sign_type HS_field_dead_spin   = configuration_e_spin[dead_spin  ].get_HS_field();
  HS_field_sign_type HS_field_living_spin = configuration_e_spin[living_spin].get_HS_field();

  std::pair<int,int>& pair_dead_spin   = full_configuration[configuration_index_dead_spin]
  .get_configuration_e_spin_indices();
  std::pair<int,int>& pair_living_spin =
  full_configuration[configuration_index_living_spin].get_configuration_e_spin_indices();

  dca::linalg::matrixop::swapRowAndCol(N , dead_spin, living_spin);
  dca::linalg::matrixop::swapRowAndCol(G0, dead_spin, living_spin);

  swap(configuration_e_spin[dead_spin], configuration_e_spin[living_spin]);

  if(HS_field_dead_spin == HS_FIELD_DN)
  pair_dead_spin.first = living_spin;
  else
  pair_dead_spin.second = living_spin;

  if(HS_field_living_spin == HS_FIELD_DN)
  pair_living_spin.first = dead_spin;
  else
  pair_living_spin.second = dead_spin;

  }

  dead_spin++;
  living_spin--;

  if(dead_spin == int(configuration_e_spin.size()) || living_spin == -1)
  break;
  else
  {
  configuration_index_dead_spin   = configuration_e_spin[dead_spin]  .get_configuration_index();
  configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
  }
  }
  }
*/

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::swap_non_changed_vertices_to_left(configuration_type& full_configuration,
                                                               std::vector<int>& source_index,
                                                               std::vector<int>& target_index,
                                                               e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  if (configuration_size == 0) {
    return;
  }

  int dead_spin = 0;
  int configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();

  int living_spin = configuration_size - 1;
  int configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();

  while (true) {
    while (dead_spin < configuration_size - 1 &&
           !(!full_configuration[configuration_index_dead_spin].is_annihilatable() &&
             !full_configuration[configuration_index_dead_spin].is_creatable())) {
      dead_spin++;
      configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
    }

    while (living_spin > 0 && !full_configuration[configuration_index_living_spin].is_creatable()) {
      living_spin--;
      configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
    }

    // if(dead_spin > living_spin)
    if (dead_spin >= living_spin)
      break;
    else {
      HS_field_sign_type HS_field_dead_spin = configuration_e_spin[dead_spin].get_HS_field();
      HS_field_sign_type HS_field_living_spin = configuration_e_spin[living_spin].get_HS_field();

      std::pair<int, int>& pair_dead_spin =
          full_configuration[configuration_index_dead_spin].get_configuration_e_spin_indices();
      std::pair<int, int>& pair_living_spin =
          full_configuration[configuration_index_living_spin].get_configuration_e_spin_indices();

      source_index.push_back(dead_spin);
      target_index.push_back(living_spin);

      // dca::linalg::matrixop::swapRowAndCol(N , dead_spin, living_spin);
      // dca::linalg::matrixop::swapRowAndCol(G0, dead_spin, living_spin);

      std::swap(configuration_e_spin[dead_spin], configuration_e_spin[living_spin]);

      if (HS_field_dead_spin == HS_FIELD_DN)
        pair_dead_spin.first = living_spin;
      else
        pair_dead_spin.second = living_spin;

      if (HS_field_living_spin == HS_FIELD_DN)
        pair_living_spin.first = dead_spin;
      else
        pair_living_spin.second = dead_spin;
    }

    dead_spin++;
    living_spin--;

    if (dead_spin == int(configuration_e_spin.size()) || living_spin == -1)
      break;
    else {
      configuration_index_dead_spin = configuration_e_spin[dead_spin].get_configuration_index();
      configuration_index_living_spin = configuration_e_spin[living_spin].get_configuration_index();
    }
  }
}

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::erase_non_creatable_and_non_annihilatable_spins(
    configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N_up,
    dca::linalg::Matrix<double, device_t>& N_dn, dca::linalg::Matrix<double, device_t>& G0_up,
    dca::linalg::Matrix<double, device_t>& G0_dn) {
  erase_non_creatable_and_non_annihilatable_spins(full_configuration, N_up, G0_up, e_UP);
  erase_non_creatable_and_non_annihilatable_spins(full_configuration, N_dn, G0_dn, e_DN);

  for (int i = 0; i < full_configuration.size();) {
    if (!full_configuration[i].is_annihilatable() && !full_configuration[i].is_creatable())
      full_configuration.remove_HS_spin(i);
    else
      i++;
  }
}

template <dca::linalg::DeviceType device_t>
template <class configuration_type>
void SHRINK_TOOLS<device_t>::erase_non_creatable_and_non_annihilatable_spins(
    configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G0, e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

  for (size_t i = 0; i < configuration_e_spin.size();) {
    int configuration_index = configuration_e_spin[i].get_configuration_index();
    if (!full_configuration[configuration_index].is_annihilatable() &&
        !full_configuration[configuration_index].is_creatable())
      configuration_e_spin.erase(configuration_e_spin.begin() + i);
    else
      i++;
  }

  //     N .size() = configuration_e_spin.size();
  //     G0.size() = configuration_e_spin.size();
  int SIZE = configuration_e_spin.size();

  N.resize(SIZE);
  G0.resize(SIZE);
}

template <dca::linalg::DeviceType device_t>
bool SHRINK_TOOLS<device_t>::test_swap_vectors(std::vector<int>& source_index,
                                               std::vector<int>& target_index, int size) {
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

  return true;
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_HPP
