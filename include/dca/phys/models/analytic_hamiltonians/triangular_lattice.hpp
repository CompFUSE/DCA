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

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TRIANGULAR_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TRIANGULAR_LATTICE_HPP

#include <complex>
#include <vector>

#include "dca/function/function.hpp"
#include "phys_library/domains/cluster/symmetries/point_groups/No_symmetry.h"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename DCA_point_group_type>
class triangular_lattice {
public:
  typedef no_symmetry<2> LDA_point_group;
  typedef DCA_point_group_type DCA_point_group;
  const static int DIMENSION = 2;
  const static int BANDS = 1;

  static double* initialize_r_DCA_basis();

  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k,
                                                   int b1, int s1, int b2, int s2);
};

template <typename DCA_point_group_type>
double* triangular_lattice<DCA_point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = cos(M_PI / 3.);
  r_DCA[1] = sin(M_PI / 3.);
  r_DCA[2] = cos(-M_PI / 3.);
  r_DCA[3] = sin(-M_PI / 3.);

  return r_DCA;
}

// template<typename DCA_point_group_type>
// double* triangular_lattice<DCA_point_group_type>::initialize_k_DCA_basis()
// {
//   static double* k_DCA = new double[4];

//   k_DCA[0] = 2*M_PI;  k_DCA[1] = 2*M_PI/sqrt(3.);
//   k_DCA[2] = 2*M_PI;  k_DCA[3] = -2*M_PI/sqrt(3.);

//   return k_DCA;
// }

template <typename DCA_point_group_type>
double* triangular_lattice<DCA_point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = cos(M_PI / 3.);
  r_LDA[1] = sin(M_PI / 3.);
  r_LDA[2] = cos(-M_PI / 3.);
  r_LDA[3] = sin(-M_PI / 3.);

  return r_LDA;
}

template <typename DCA_point_group_type>
std::vector<int> triangular_lattice<DCA_point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;

  return flavors;
}

template <typename DCA_point_group_type>
std::vector<std::vector<double>> triangular_lattice<DCA_point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename DCA_point_group_type>
template <class domain, class parameters_type>
void triangular_lattice<DCA_point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_interaction, parameters_type& parameters) {
  double U = parameters.get_U();

  H_interaction(0, 0, 0) = 0;
  H_interaction(0, 1, 0) = U;
  H_interaction(1, 0, 0) = U;
  H_interaction(1, 1, 0) = 0;
}

template <typename DCA_point_group_type>
template <class domain>
void triangular_lattice<DCA_point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries(0, 0) = 0;
  H_symmetries(0, 1) = -1;
  H_symmetries(1, 0) = -1;
  H_symmetries(1, 1) = 0;
}

template <typename DCA_point_group_type>
template <class parameters_type>
std::complex<double> triangular_lattice<DCA_point_group_type>::get_LDA_Hamiltonians(
    parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2) {
  assert(k.size() == 2);

  std::complex<double> H_LDA = 0.;

  double t = parameters.get_t();

  if ((b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1))) {
    H_LDA = 0;

    H_LDA += -2. * t * cos(k[0]);
    H_LDA += -4. * t * cos(sqrt(3.) * k[1] / 2.) * cos(k[0] / 2.);
  }

  return H_LDA;
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TRIANGULAR_LATTICE_HPP
