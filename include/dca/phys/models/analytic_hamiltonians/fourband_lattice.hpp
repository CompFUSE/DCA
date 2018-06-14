// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// 4-band lattice.
//
// TODO: - Replace get_LDA_Hamiltonians with intialize_H_0 function (see e.g. square_lattice.hpp).
//       - Use correct index of origin in initialize_H_interaction (see e.g. square_lattice.hpp).

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FOURBAND_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FOURBAND_LATTICE_HPP

#include <cmath>
#include <complex>
#include <utility>
#include <vector>

#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class fourband_lattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 4;

public:
  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <class domain, class parameters_type>
  static void initialize_H_LDA(func::function<std::complex<double>, domain>& H_LDA,
                               parameters_type& parameters);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters,
                                                   std::vector<double> /*k*/, int b1, int s1,
                                                   int b2, int s2);
};

template <typename point_group_type>
double* fourband_lattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* fourband_lattice<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* fourband_lattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* fourband_lattice<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> fourband_lattice<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> fourband_lattice<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> fourband_lattice<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <class domain, class parameters_type>
void fourband_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_interaction, parameters_type& parameters) {
  H_interaction = 0.;

  double U0 = parameters.get_U0();
  double U1 = parameters.get_U1();
  double V = parameters.get_V();
  double V_prime = parameters.get_V_prime();

  H_interaction(0, 0, 0, 1, 0) = U0;
  H_interaction(0, 1, 0, 0, 0) = U0;

  H_interaction(1, 0, 1, 1, 0) = U1;
  H_interaction(1, 1, 1, 0, 0) = U1;

  H_interaction(1, 0, 0, 1, 0) = V;
  H_interaction(0, 1, 1, 0, 0) = V;
  H_interaction(0, 0, 1, 1, 0) = V;
  H_interaction(1, 1, 0, 0, 0) = V;

  H_interaction(1, 0, 0, 0, 0) = V_prime;
  H_interaction(0, 0, 1, 0, 0) = V_prime;
  H_interaction(0, 1, 1, 1, 0) = V_prime;
  H_interaction(1, 1, 0, 1, 0) = V_prime;
}

template <typename point_group_type>
template <class domain>
void fourband_lattice<point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;

  H_symmetries(2, 0, 2, 0) = 2;
  H_symmetries(2, 1, 2, 1) = 2;

  H_symmetries(3, 0, 3, 0) = 3;
  H_symmetries(3, 1, 3, 1) = 3;
}

template <typename point_group_type>
template <class domain, class parameters_type>
void fourband_lattice<point_group_type>::initialize_H_LDA(
    func::function<std::complex<double>, domain>& H_LDA, parameters_type& parameters) {
  typedef typename parameters_type::b b;
  typedef typename parameters_type::s s;

  typedef typename parameters_type::k_LDA k_LDA;

  std::vector<double> k;

  for (int k_ind = 0; k_ind < k_LDA::dmn_size(); k_ind++) {
    k = k_LDA::parameter_type::get_elements()[k_ind];

    for (int b_ind1 = 0; b_ind1 < b::dmn_size(); b_ind1++)
      for (int s_ind1 = 0; s_ind1 < s::dmn_size(); s_ind1++)
        for (int b_ind2 = 0; b_ind2 < b::dmn_size(); b_ind2++)
          for (int s_ind2 = 0; s_ind2 < s::dmn_size(); s_ind2++)
            H_LDA(b_ind1, s_ind1, b_ind2, s_ind2, k_ind) =
                get_LDA_Hamiltonians(parameters, k, b_ind1, s_ind1, b_ind2, s_ind2);
  }
}

template <typename point_group_type>
template <class parameters_type>
std::complex<double> fourband_lattice<point_group_type>::get_LDA_Hamiltonians(
    parameters_type& parameters, std::vector<double> /*k*/, int b1, int s1, int b2, int s2) {
  const static std::complex<double> I(0, 1);

  double ei0 = parameters.get_ei0();
  double eb0 = parameters.get_eb0();
  double t0 = parameters.get_t0();

  double ei1 = parameters.get_ei1();
  double eb1 = parameters.get_eb1();
  double t1 = parameters.get_t1();

  std::complex<double> H_LDA = 0.;

  if (s1 == s2) {
    if (b1 == b2) {
      if (b1 == 0)
        H_LDA += ei0;
      else if (b1 == 1)
        H_LDA += ei1;
      else if (b1 == 2)
        H_LDA += eb0;
      else
        H_LDA += eb1;
    }
    else if (abs(b1 - b2) == 2) {
      if (b1 % 2 == 0)
        H_LDA += t0;
      else if (b1 % 2 == 1)
        H_LDA += t1;
    }
  }

  return H_LDA;
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FOURBAND_LATTICE_HPP
