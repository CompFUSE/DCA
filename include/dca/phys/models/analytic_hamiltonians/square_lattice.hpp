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

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP

#include <complex>
#include <utility>
#include <vector>

#include "dca/function/function.hpp"
#include "dca/util/type_list.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class square_lattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef point_group_type DCA_point_group;

  const static int DIMENSION = 2;
  const static int BANDS = 1;

  static double* initialize_r_DCA_basis();

  static double* initialize_r_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  template <class domain, class parameters_type>
  static void initialize_H_interaction(func::function<double, domain>& H_interaction,
                                       parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  template <class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k,
                                                   int b1, int s1, int b2, int s2);
};

template <typename point_group_type>
double* square_lattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.;
  r_DCA[1] = 0.;
  r_DCA[2] = 0.;
  r_DCA[3] = 1.;

  return r_DCA;
}

template <typename point_group_type>
double* square_lattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
std::vector<int> square_lattice<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> square_lattice<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> square_lattice<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

// TODO: Add non-local interaction of same spins.
//       Use V instead of U_prime?
template <typename point_group_type>
template <class domain, class parameters_type>
void square_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, domain>& H_interaction, parameters_type& parameters) {
  H_interaction = 0.;

  // actually the same as DCA_r_cluster_type (see typedifinitions.h).
  typedef
      typename dca::util::TypeAt<0, typename domain::template domain_typelist<2>>::type DCA_r_cluster_t;

  int DIMENSION = DCA_r_cluster_t::DIMENSION;
  assert(DIMENSION == 2);

  int origin = DCA_r_cluster_t::origin_index();

  std::vector<typename DCA_r_cluster_t::element_type>& basis = DCA_r_cluster_t::get_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& super_basis =
      DCA_r_cluster_t::get_super_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& elements = DCA_r_cluster_t::get_elements();

  std::vector<int> nn_index(DIMENSION);  // Indices of nearest neighbours w.r.t. origin.
  for (int d = 0; d < DIMENSION; ++d) {
    std::vector<double> basis_vec =
        domains::cluster_operations::translate_inside_cluster(basis[d], super_basis);
    nn_index[d] = domains::cluster_operations::index(basis_vec, elements, domains::BRILLOUIN_ZONE);
  }

  // Nearest-neighbor opposite spin interaction
  double V = parameters.get_V();
  H_interaction(0, 1, nn_index[0]) = V;
  H_interaction(1, 0, nn_index[0]) = V;
  H_interaction(0, 1, nn_index[1]) = V;
  H_interaction(1, 0, nn_index[1]) = V;

  // Nearest-neighbor same spin interaction
  double V_prime = parameters.get_V_prime();
  H_interaction(0, 0, nn_index[0]) = V_prime;
  H_interaction(1, 1, nn_index[0]) = V_prime;
  H_interaction(0, 0, nn_index[1]) = V_prime;
  H_interaction(1, 1, nn_index[1]) = V_prime;

  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  double U = parameters.get_U();
  H_interaction(0, 1, origin) = U;
  H_interaction(1, 0, origin) = U;
}

template <typename point_group_type>
template <class domain>
void square_lattice<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries(0, 0) = 0;
  H_symmetries(0, 1) = -1;
  H_symmetries(1, 0) = -1;
  H_symmetries(1, 1) = 0;
}

template <typename point_group_type>
template <class parameters_type>
std::complex<double> square_lattice<point_group_type>::get_LDA_Hamiltonians(
    parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2) {
  std::complex<double> H_LDA = 0.;

  double t = parameters.get_t();
  double t_prime = parameters.get_t_prime();

  if ((b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2. * t * (cos(k[0]) + cos(k[1])) - 4. * t_prime * cos(k[0]) * cos(k[1]);

  return H_LDA;
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_SQUARE_LATTICE_HPP
