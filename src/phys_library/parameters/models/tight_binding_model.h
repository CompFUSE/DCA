// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_PARAMETERS_MODELS_TIGHT_BINDING_MODEL_H
#define PHYS_LIBRARY_PARAMETERS_MODELS_TIGHT_BINDING_MODEL_H

#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

template <typename lattice>
class tight_binding_model {
public:
  using lattice_type = lattice;
  using k_LDA =
      dmn_0<cluster_domain<double, lattice_type::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, PARALLELLEPIPEDUM>>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;

public:
  // taken care off via parameters !
  static int& get_DCA_size();
  static int& get_LDA_size();

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  // taken care of via lattice-type
  typedef typename lattice_type::LDA_point_group LDA_point_group;
  typedef typename lattice_type::DCA_point_group DCA_point_group;

  static const int DIMENSION = lattice_type::DIMENSION;
  static const int BANDS = lattice_type::BANDS;

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_interactions_symmetries);

  template <class domain, class parameters_type>
  static void initialize_H_LDA(FUNC_LIB::function<std::complex<double>, domain>& H_LDA,
                               parameters_type& parameters);

  template <class parameters_type>
  static void initialize(parameters_type& parameters);

private:
  template <class parameters_type>
  static void initialize_Bett_cluster(parameters_type& parameters);

  template <class parameters_type>
  static void initialize_default(parameters_type& parameters);
};

template <typename lattice_type>
int& tight_binding_model<lattice_type>::get_DCA_size() {
  static int DCA_size = -1;
  return DCA_size;
}

template <typename lattice_type>
int& tight_binding_model<lattice_type>::get_LDA_size() {
  static int LDA_size = -1;
  return LDA_size;
}

template <typename lattice_type>
std::vector<int>& tight_binding_model<lattice_type>::DCA_grid_size() {
  static std::vector<int> v(0);
  return v;
}

template <typename lattice_type>
std::vector<int>& tight_binding_model<lattice_type>::LDA_grid_size() {
  static std::vector<int> v(0);
  return v;
}

template <typename lattice_type>
double* tight_binding_model<lattice_type>::get_r_DCA_basis() {
  static double* r_DCA = lattice_type::initialize_r_DCA_basis();
  return r_DCA;
}

template <typename lattice_type>
double* tight_binding_model<lattice_type>::get_r_LDA_basis() {
  static double* r_LDA = lattice_type::initialize_r_LDA_basis();
  return r_LDA;
}

template <typename lattice_type>
std::vector<int> tight_binding_model<lattice_type>::get_flavors() {
  return lattice_type::get_flavors();
}

template <typename lattice_type>
std::vector<std::vector<double>> tight_binding_model<lattice_type>::get_a_vectors() {
  return lattice_type::get_a_vectors();
}

template <typename lattice_type>
template <class domain, class parameters_type>
void tight_binding_model<lattice_type>::initialize_H_interaction(
    FUNC_LIB::function<double, domain>& H_interaction, parameters_type& parameters) {
  lattice_type::initialize_H_interaction(H_interaction, parameters);
}

template <typename lattice_type>
template <class domain>
void tight_binding_model<lattice_type>::initialize_H_symmetries(
    FUNC_LIB::function<int, domain>& H_symmetry) {
  lattice_type::initialize_H_symmetry(H_symmetry);
}

template <typename lattice_type>
template <class domain, class parameters_type>
void tight_binding_model<lattice_type>::initialize_H_LDA(
    FUNC_LIB::function<std::complex<double>, domain>& H_LDA, parameters_type& parameters) {
  std::vector<double> k;

  for (int k_ind = 0; k_ind < k_LDA::dmn_size(); k_ind++) {
    k = k_LDA::get_elements()[k_ind];

    for (int b_ind1 = 0; b_ind1 < b::dmn_size(); b_ind1++)
      for (int s_ind1 = 0; s_ind1 < s::dmn_size(); s_ind1++)
        for (int b_ind2 = 0; b_ind2 < b::dmn_size(); b_ind2++)
          for (int s_ind2 = 0; s_ind2 < s::dmn_size(); s_ind2++)
            H_LDA(b_ind1, s_ind1, b_ind2, s_ind2, k_ind) =
                lattice_type::get_LDA_Hamiltonians(parameters, k, b_ind1, s_ind1, b_ind2, s_ind2);
  }
}

template <typename lattice_type>
template <class parameters_type>
void tight_binding_model<lattice_type>::initialize(parameters_type& /*parameters*/) {}

#endif  // PHYS_LIBRARY_PARAMETERS_MODELS_TIGHT_BINDING_MODEL_H
