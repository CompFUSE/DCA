// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Tight binding model.

#ifndef DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP
#define DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename Lattice>
class TightBindingModel {
public:
  using lattice_type = Lattice;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

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
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetries(func::function<int, domain>& H_interactions_symmetries);

  template <class parameters_type>
  static void initialize(parameters_type& /*parameters*/) {}

  template <typename ParametersType, typename FunctionType>
  static void initialize_H_0(const ParametersType& parameters, FunctionType& H_0) {
    Lattice::initialize_H_0(parameters, H_0);
  }

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
void TightBindingModel<Lattice>::initialize_H_interaction(func::function<double, domain>& H_interaction,
                                                          parameters_type& parameters) {
  Lattice::initialize_H_interaction(H_interaction, parameters);
}

template <typename Lattice>
template <class domain>
void TightBindingModel<Lattice>::initialize_H_symmetries(func::function<int, domain>& H_symmetry) {
  Lattice::initialize_H_symmetry(H_symmetry);
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP
