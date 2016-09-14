// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP
#define DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP

#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename Lattice>
class TightBindingModel {
public:
  using lattice_type = Lattice;
  using k_LDA =
      dmn_0<cluster_domain<double, Lattice::DIMENSION, LATTICE_SP, MOMENTUM_SPACE, PARALLELLEPIPEDUM>>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;

public:
  // taken care off via parameters !
  static int& get_DCA_size();
  static int& get_LDA_size();

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  // taken care of via lattice-type
  typedef typename Lattice::LDA_point_group LDA_point_group;
  typedef typename Lattice::DCA_point_group DCA_point_group;

  static const int DIMENSION = Lattice::DIMENSION;
  static const int BANDS = Lattice::BANDS;

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

template <typename Lattice>
int& TightBindingModel<Lattice>::get_DCA_size() {
  static int DCA_size = -1;
  return DCA_size;
}

template <typename Lattice>
int& TightBindingModel<Lattice>::get_LDA_size() {
  static int LDA_size = -1;
  return LDA_size;
}

template <typename Lattice>
std::vector<int>& TightBindingModel<Lattice>::DCA_grid_size() {
  static std::vector<int> v(0);
  return v;
}

template <typename Lattice>
std::vector<int>& TightBindingModel<Lattice>::LDA_grid_size() {
  static std::vector<int> v(0);
  return v;
}

template <typename Lattice>
double* TightBindingModel<Lattice>::get_r_DCA_basis() {
  static double* r_DCA = Lattice::initialize_r_DCA_basis();
  return r_DCA;
}

template <typename Lattice>
double* TightBindingModel<Lattice>::get_r_LDA_basis() {
  static double* r_LDA = Lattice::initialize_r_LDA_basis();
  return r_LDA;
}

template <typename Lattice>
std::vector<int> TightBindingModel<Lattice>::get_flavors() {
  return Lattice::get_flavors();
}

template <typename Lattice>
std::vector<std::vector<double>> TightBindingModel<Lattice>::get_a_vectors() {
  return Lattice::get_a_vectors();
}

template <typename Lattice>
template <class domain, class parameters_type>
void TightBindingModel<Lattice>::initialize_H_interaction(
    FUNC_LIB::function<double, domain>& H_interaction, parameters_type& parameters) {
  Lattice::initialize_H_interaction(H_interaction, parameters);
}

template <typename Lattice>
template <class domain>
void TightBindingModel<Lattice>::initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_symmetry) {
  Lattice::initialize_H_symmetry(H_symmetry);
}

template <typename Lattice>
template <class domain, class parameters_type>
void TightBindingModel<Lattice>::initialize_H_LDA(
    FUNC_LIB::function<std::complex<double>, domain>& H_LDA, parameters_type& parameters) {
  std::vector<double> k;

  for (int k_ind = 0; k_ind < k_LDA::dmn_size(); k_ind++) {
    k = k_LDA::get_elements()[k_ind];

    for (int b_ind1 = 0; b_ind1 < b::dmn_size(); b_ind1++)
      for (int s_ind1 = 0; s_ind1 < s::dmn_size(); s_ind1++)
        for (int b_ind2 = 0; b_ind2 < b::dmn_size(); b_ind2++)
          for (int s_ind2 = 0; s_ind2 < s::dmn_size(); s_ind2++)
            H_LDA(b_ind1, s_ind1, b_ind2, s_ind2, k_ind) =
                Lattice::get_LDA_Hamiltonians(parameters, k, b_ind1, s_ind1, b_ind2, s_ind2);
  }
}

template <typename Lattice>
template <class parameters_type>
void TightBindingModel<Lattice>::initialize(parameters_type& /*parameters*/) {}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP
