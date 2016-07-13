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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_SP_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_SP_H

#include <complex>
#include <iostream>
#include <sstream>
#include <vector>

#include "dca/math/geometry/gaussian_quadrature/gaussian_quadrature_domain.hpp"
#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "math_library/geometry_library/tetrahedron_mesh/tetrahedron_mesh.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_interpolation_matrices.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_routines.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/tetrahedron_integration.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/tetrahedron_routines_harmonic_function.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"

namespace DCA {

template <typename parameters_type, typename K_dmn>
class coarsegraining_sp : public coarsegraining_routines<parameters_type, K_dmn>,
                          public tetrahedron_integration<parameters_type, K_dmn> {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using k_cluster_type = typename K_dmn::parameter_type;
  using host_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                             LATTICE_SP, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_HOST = dmn_0<host_k_cluster_type>;

#ifdef SINGLE_PRECISION_COARSEGRAINING
  using scalar_type = float;
#else
  using scalar_type = double;
#endif  // SINGLE_PRECISION_COARSEGRAINING
  using complex_type = std::complex<scalar_type>;

  using tetrahedron_dmn = dmn_0<math_algorithms::tetrahedron_mesh<K_dmn>>;
  using quadrature_dmn = dca::math::gaussquad::gaussian_quadrature_domain<tetrahedron_dmn>;

  using q_dmn = dmn_0<coarsegraining_domain<K_dmn, K>>;
  using q_0_dmn = dmn_0<coarsegraining_domain<K_dmn, ORIGIN>>;
  using tet_dmn = dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K>>;
  using tet_0_dmn = dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>>;

  using t = dmn_0<time_domain>;
  using w = dmn_0<frequency_domain>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  using nu_nu_q = dmn_variadic<nu, nu, q_dmn>;
  using nu_nu_tet = dmn_variadic<nu, nu, tet_dmn>;

  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

public:
  coarsegraining_sp(parameters_type& parameters_ref);

  void initialize();

  void reset_fine_q_mesh(int recursion, int rule, int period);

  template <typename other_scalar_type, typename r_dmn>
  void compute_phi_r(FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r);

  template <typename other_scalar_type, typename r_dmn>
  void print_phi_r_to_shell(FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r);

  /*****************************************
   ***                                   ***
   ***         Routines for DCA          ***
   ***                                   ***
   *****************************************/

  template <typename other_scalar_type>
  void compute_S_K_w(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_k_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w);

  template <typename other_scalar_type, typename k_dmn>
  void compute_G_K_w(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& G_K_w);

  template <typename other_scalar_type, typename k_dmn, typename w_dmn>
  void compute_G_K_w_with_TIM(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>>& S_K_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>>& G_K_w);

  /*****************************************
   ***                                   ***
   ***         Routines for DCA+         ***
   ***                                   ***
   *****************************************/

  template <typename other_scalar_type, typename k_dmn>
  void compute_S_K_w(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>>& S_k_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w);

  template <typename other_scalar_type, typename k_dmn>
  void compute_G_K_w(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>>& S_k_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& G_K_w);

  template <typename other_scalar_type, typename k_dmn, typename w_dmn>
  void compute_G_K_w_with_TIM(
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w_dmn>>& S_K_w,
      FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>>& G_K_w);

  template <typename other_scalar_type, typename k_dmn>
  void compute_G0_K_t(FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
                      FUNC_LIB::function<other_scalar_type, dmn_4<nu, nu, K_dmn, t>>& G0_k_w);

private:
  void update_shell(int i, int N);

  template <typename other_scalar_type, typename r_dmn>
  void plot_phi_r(FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  bool initialized;

  // tetrahedron q-points
  FUNC_LIB::function<scalar_type, tet_dmn> w_tet;

  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_tet> I_tet;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_tet> H_tet;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_tet> S_tet;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_tet> A_tet;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_tet> G_tet;

  // gaussian q-points
  FUNC_LIB::function<scalar_type, q_dmn> w_q;

  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> I_q;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> H_q;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> S_q;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> A_q;
  FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> G_q;
};

template <typename parameters_type, typename K_dmn>
coarsegraining_sp<parameters_type, K_dmn>::coarsegraining_sp(parameters_type& parameters_ref)
    : coarsegraining_routines<parameters_type, K_dmn>(parameters_ref),
      tetrahedron_integration<parameters_type, K_dmn>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      initialized(false),

      // tetrahedron q-points
      w_tet("w_tet"),

      I_tet("I_tet"),
      H_tet("H_tet"),
      S_tet("S_tet"),
      A_tet("A_tet"),
      G_tet("G_tet"),

      // gaussian q-points
      w_q("w_q"),

      I_q("I_q"),
      H_q("H_q"),
      S_q("S_q"),
      A_q("A_q"),
      G_q("G_q") {
  initialize();
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::update_shell(int i, int N) {
  int tmp = i;

  if (concurrency.id() == concurrency.first() && N > 10 && (tmp % (N / 10)) == 0) {
    std::stringstream ss;

    ss << std::scientific;
    ss.precision(6);

    ss << "\t\t\t" << double(i) / double(N) * 100. << " % completed \t ";
    ss << dca::util::print_time();
    ss << "\n";

    std::cout << ss.str();
  }
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::initialize() {
  {
    this->compute_tetrahedron_mesh(parameters.get_k_mesh_refinement(),
                                   parameters.get_number_of_periods());

    this->compute_gaussian_mesh(parameters.get_k_mesh_refinement(),
                                parameters.get_gaussian_quadrature_rule(),
                                parameters.get_number_of_periods());
  }

  {
    w_q.reset();
    w_tet.reset();

    for (int l = 0; l < w_q.size(); l++)
      w_q(l) = q_dmn::parameter_type::get_weights()[l];

    for (int l = 0; l < w_tet.size(); l++)
      w_tet(l) = tet_dmn::parameter_type::get_weights()[l];
  }

  {
    I_tet.reset();
    H_tet.reset();
    S_tet.reset();
    A_tet.reset();
    G_tet.reset();

    I_q.reset();
    H_q.reset();
    S_q.reset();
    A_q.reset();
    G_q.reset();
  }

  interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename r_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::plot_phi_r(
    FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r) {
  std::vector<double> x;
  std::vector<double> y;

  std::vector<std::vector<double>> super_basis = r_dmn::parameter_type::get_super_basis_vectors();
  for (int r_ind = 0; r_ind < r_dmn::dmn_size(); r_ind++) {
    std::vector<double> r_vec = r_dmn::get_elements()[r_ind];
    std::vector<std::vector<double>> r_vecs =
        cluster_operations::equivalent_vectors(r_vec, super_basis);

    x.push_back(std::sqrt(r_vecs[0][0] * r_vecs[0][0] + r_vecs[0][1] * r_vecs[0][1]));
    y.push_back((phi_r(r_ind)));
  }

  SHOW::plot_points(x, y);
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename r_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_phi_r(
    FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r) {
  math_algorithms::tetrahedron_mesh<k_cluster_type> mesh(parameters.get_k_mesh_refinement());

  quadrature_dmn::translate_according_to_period(parameters.get_number_of_periods(), mesh);

  std::vector<math_algorithms::tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  {
    phi_r = 0.;

    r_dmn r_domain;
    std::pair<int, int> bounds = concurrency.get_bounds(r_domain);

    std::vector<std::vector<double>> super_basis = r_dmn::parameter_type::get_super_basis_vectors();

    for (int l = bounds.first; l < bounds.second; l++) {
      std::vector<double> r_vec = r_dmn::get_elements()[l];
      std::vector<std::vector<double>> r_vecs =
          cluster_operations::equivalent_vectors(r_vec, super_basis);
      for (int r_ind = 0; r_ind < r_vecs.size(); r_ind++)
        for (int tet_ind = 0; tet_ind < tetrahedra.size(); tet_ind++)
          phi_r(l) += std::real(tetrahedron_routines_harmonic_function::execute(
                          r_vecs[0], tetrahedra[tet_ind])) /
                      r_vecs.size();
      //        }
    }

    concurrency.sum(phi_r);

    {
      scalar_type V_K = 0;
      for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      phi_r /= V_K;
    }
  }
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename r_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::print_phi_r_to_shell(
    FUNC_LIB::function<other_scalar_type, r_dmn>& phi_r) {
  std::cout << "\n" << __FUNCTION__ << std::endl;

  math_algorithms::tetrahedron_mesh<k_cluster_type> mesh(parameters.get_k_mesh_refinement());

  quadrature_dmn::translate_according_to_period(parameters.get_number_of_periods(), mesh);

  std::vector<math_algorithms::tetrahedron<DIMENSION>>& tetrahedra = mesh.get_tetrahedra();

  phi_r = 0.;

  r_dmn r_domain;
  std::pair<int, int> bounds = concurrency.get_bounds(r_domain);

  for (int l = bounds.first; l < bounds.second; l++) {
    std::vector<double> r_vec = r_dmn::get_elements()[l];

    for (int tet_ind = 0; tet_ind < tetrahedra.size(); tet_ind++)
      phi_r(l) +=
          std::real(tetrahedron_routines_harmonic_function::execute(r_vec, tetrahedra[tet_ind]));
  }

  concurrency.sum(phi_r);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    phi_r /= V_K;
  }

  std::cout << "r_ind\tx\t\ty\t\tr\t\tphi(r)" << std::endl;

  for (int r_ind = 0; r_ind < r_dmn::dmn_size(); r_ind++) {
    std::vector<double> r_vec = r_dmn::get_elements()[r_ind];

    std::cout << r_ind << "\t" << r_vec[0] << "\t" << r_vec[1] << "\t"
              << std::sqrt(r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1]) << "\t" << phi_r(r_ind)
              << std::endl;
  }
}

/*****************************************
 ***                                   ***
 ***         Routines for DCA          ***
 ***                                   ***
 *****************************************/

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type>
void coarsegraining_sp<parameters_type, K_dmn>::compute_S_K_w(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_k_w,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w) {
  for (int l = 0; l < S_k_w.size(); l++)
    S_K_w(l) = S_k_w(l);
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename k_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& G_K_w) {
  // cout << "\n\n\n\t" << __FUNCTION__ << "\n\n\n";

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> H_k("H_k");
  FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w>> S_K("S_K");

  for (int l = 0; l < S_K_w.size(); l++)
    S_K(l) = S_K_w(l);

  for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        H_k(i, j, k_ind) = H_0(i, j, k_ind);

  G_K_w = 0.;

  dmn_2<K_dmn, w> K_wm_dmn;
  std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

  int coor[2];
  for (int l = bounds.first; l < bounds.second; l++) {
    K_wm_dmn.linind_2_subind(l, coor);

    this->compute_G_q_w(coor[0], coor[1], H_k, S_K, I_q, H_q, S_q, G_q);

    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_K_w(i, j, coor[0], coor[1]) += G_q(i, j, q_ind) * w_q(q_ind);
  }

  concurrency.sum(G_K_w);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    G_K_w /= V_K;
  }
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename k_dmn, typename w_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w_with_TIM(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>>& S_K_w,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>>& G_K_w) {
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> H_k("H_k");
  FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn>> S_K("S_K");

  for (int l = 0; l < S_K_w.size(); l++)
    S_K(l) = S_K_w(l);

  for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        H_k(i, j, k_ind) = H_0(i, j, k_ind);

  G_K_w = 0.;

  dmn_2<K_dmn, w_dmn> K_wm_dmn;
  std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

  int coor[2];
  for (int l = bounds.first; l < bounds.second; l++) {
    K_wm_dmn.linind_2_subind(l, coor);

    this->compute_G_q_w(coor[0], coor[1], H_k, S_K, I_tet, H_tet, S_tet, G_tet);

    {
      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu>> G_int;

      this->tetrahedron_integration_mt(w_tet, G_tet, G_int);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_K_w(i, j, coor[0], coor[1]) += G_int(i, j);
    }
  }

  concurrency.sum(G_K_w);

  {
    scalar_type V_K = 0;
    for (int tet_ind = 0; tet_ind < tet_dmn::dmn_size(); tet_ind++)
      V_K += w_tet(tet_ind);

    G_K_w /= V_K;
  }
}

/*****************************************
 ***                                   ***
 ***         Routines for DCA+         ***
 ***                                   ***
 *****************************************/

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename k_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_S_K_w(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>>& S_k_w,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& S_K_w) {
  S_K_w = 0.;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> A_k("A_k");
  FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>> A_k_w("A_k_w");

  transform_to_alpha::forward(1., S_k_w, A_k_w);

  dmn_2<K_dmn, w> K_wm_dmn;
  std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

  int coor[2];
  for (int l = bounds.first; l < bounds.second; l++) {
    K_wm_dmn.linind_2_subind(l, coor);

    for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          A_k(i, j, k_ind) = A_k_w(i, j, k_ind, coor[1]);

    this->compute_S_q_from_A_k(coor[0], coor[1], A_k, A_q, S_q);

    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_K_w(i, j, coor[0], coor[1]) += S_q(i, j, q_ind) * w_q(q_ind);
  }

  concurrency.sum(S_K_w);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    S_K_w /= V_K;
  }
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename k_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>>& S_k_w,
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w>>& G_K_w) {
  G_K_w = 0.;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> H_k("H_k");
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> A_k("A_k");

  FUNC_LIB::function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w>> A_k_w("A_k_w");

  transform_to_alpha::forward(1., S_k_w, A_k_w);

  for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        H_k(i, j, k_ind) = H_0(i, j, k_ind);

  dmn_2<K_dmn, w> K_wm_dmn;
  std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

  int coor[2];
  for (int l = bounds.first; l < bounds.second; l++) {
    K_wm_dmn.linind_2_subind(l, coor);

    for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          A_k(i, j, k_ind) = A_k_w(i, j, k_ind, coor[1]);

    this->compute_G_q_w(coor[0], coor[1], H_k, A_k, I_q, H_q, A_q, S_q, G_q);

    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_K_w(i, j, coor[0], coor[1]) += G_q(i, j, q_ind) * w_q(q_ind);
  }

  concurrency.sum(G_K_w);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    G_K_w /= V_K;
  }
}

template <typename parameters_type, typename K_dmn>
template <typename other_scalar_type, typename k_dmn>
void coarsegraining_sp<parameters_type, K_dmn>::compute_G0_K_t(
    FUNC_LIB::function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn>>& H_0,
    FUNC_LIB::function<other_scalar_type, dmn_4<nu, nu, K_dmn, t>>& G_K_t) {
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn>> H_k("H_k");

  for (int k_ind = 0; k_ind < k_dmn::dmn_size(); k_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        H_k(i, j, k_ind) = H_0(i, j, k_ind);

  G_K_t = 0.;

  dmn_2<K_dmn, t> K_t_dmn;
  std::pair<int, int> bounds = concurrency.get_bounds(K_t_dmn);

  int coor[2];
  for (int l = bounds.first; l < bounds.second; l++) {
    K_t_dmn.linind_2_subind(l, coor);

    this->compute_G_q_t(coor[0], coor[1], H_k, I_q, H_q, G_q);

    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_K_t(i, j, coor[0], coor[1]) += std::real(G_q(i, j, q_ind)) * w_q(q_ind);
  }

  concurrency.sum(G_K_t);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    G_K_t /= V_K;
  }
}
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_SP_H
