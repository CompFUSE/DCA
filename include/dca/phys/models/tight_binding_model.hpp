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
#include <type_traits>

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
  using model_type = TightBindingModel<Lattice>;
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

  static const double* get_r_DCA_basis();
  static const double* get_k_DCA_basis();

  static const double* get_r_LDA_basis();
  static const double* get_k_LDA_basis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  template <class domain, class parameters_type>
  static void initializeHInteraction(func::function<typename parameters_type::Real, domain>& H_interaction,
                                     parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetries(func::function<int, domain>& H_interactions_symmetries);

  template <class parameters_type>
  static void initialize(parameters_type& /*parameters*/) {}

  template <typename ParametersType, typename FunctionType>
  static void initializeH0(const ParametersType& parameters, FunctionType& H_0) {
    Lattice::initializeH0(parameters, H_0);
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
const double* TightBindingModel<Lattice>::get_r_DCA_basis() {
  if constexpr (std::is_same_v<decltype(Lattice::initializeRDCABasis()), const double *>) {
    static const double* r_DCA = Lattice::initializeRDCABasis();
    return r_DCA;
  }
  else {
    static std::array<double, 9> r_DCA = Lattice::initializeRDCABasis();
    return r_DCA.data();
  }
}

template <typename Lattice>
const double* TightBindingModel<Lattice>::get_r_LDA_basis() {
  if constexpr (std::is_same_v<decltype(Lattice::initializeRLDABasis()), const double *>) {
    static const double* r_LDA = Lattice::initializeRLDABasis();
    return r_LDA;
  }
  else {
    static std::array<double, 9> r_LDA = Lattice::initializeRLDABasis();
    return r_LDA.data();
  }
}

template <typename Lattice>
std::vector<int> TightBindingModel<Lattice>::flavors() {
  return Lattice::flavors();
}

template <typename Lattice>
std::vector<std::vector<double>> TightBindingModel<Lattice>::aVectors() {
  return Lattice::aVectors();
}

template <typename Lattice>
template <class domain, class parameters_type>
void TightBindingModel<Lattice>::initializeHInteraction(func::function<typename parameters_type::Real, domain>& H_interaction,
                                                        parameters_type& parameters) {
  Lattice::initializeHInteraction(H_interaction, parameters);
}

template <typename Lattice>
template <class domain>
void TightBindingModel<Lattice>::initialize_H_symmetries(func::function<int, domain>& H_symmetry) {
  Lattice::initializeHSymmetry(H_symmetry);
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_TIGHT_BINDING_MODEL_HPP
