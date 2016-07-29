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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_ROUTINES_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_ROUTINES_H

#include <complex>

#include "dca/math/geometry/gaussian_quadrature/gaussian_quadrature_domain.hpp"

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"

#include "math_library/geometry_library/tetrahedron_mesh/tetrahedron_mesh.h"

#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_domain.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/quadrature_integration.h"
#include "phys_library/DCA+_step/lattice_mapping/interpolation/interpolation_matrices.h"
#include "phys_library/DCA+_step/lattice_mapping/interpolation/transform_to_alpha.hpp"

#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"

namespace DCA {

template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
struct coarsegraining_functions {
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  int K_ind;
  int w_ind;

  scalar_type mu;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>* H_k;
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>* A_k;

  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>> I_q;
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>> H_q;
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>> A_q;
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>> S_q;
  FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>> G_q;
};

template <typename parameters_type, typename K_dmn>
class coarsegraining_routines {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using k_cluster_type = typename K_dmn::parameter_type;

  using tetrahedron_dmn = dmn_0<math_algorithms::tetrahedron_mesh<K_dmn>>;
  using quadrature_dmn = dca::math::gaussquad::gaussian_quadrature_domain<tetrahedron_dmn>;

  using t = dmn_0<time_domain>;
  using w = dmn_0<frequency_domain>;
  using w_REAL = dmn_0<frequency_domain_real_axis>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = dmn_variadic<nu, nu>;

  using DCA_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_DCA = dmn_0<DCA_k_cluster_type>;

  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

public:
  coarsegraining_routines(parameters_type& parameters_ref);

protected:
  void compute_tetrahedron_mesh(int k_mesh_refinement, int number_of_periods);

  void compute_gaussian_mesh(int k_mesh_refinement, int gaussian_quadrature_rule,
                             int number_of_periods);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void wannier_interpolation(
      int K_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& f_k,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& f_q);

  template <typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
  void compute_I_q(tmp_scalar_type value,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_H_q(int K_ind,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q);

  /***********************************************************************************
   ***
   ***         Routines for DCA
   ***
   ***********************************************************************************/

  template <typename scalar_type, typename q_dmn_t>
  void compute_S_q(int K_ind, int w_ind,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w>>& S_K_w,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(int K_ind, int w_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w>>& S_K,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type, typename q_dmn_t>
  void compute_S_q(int K_ind, int w_ind,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w_REAL>>& S_K_w,
                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(int K_ind, int w_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w_REAL>>& S_K,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_t(int K_ind, int t_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  /***********************************************************************************
   ***                                                                             ***
   ***         Routines for DCA+
   ***                                                                             ***
   ***********************************************************************************/

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_S_q_from_A_k(
      int K_ind, int w_ind,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& A_k,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& A_q,
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(int K_ind, int w_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& A_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& A_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q);

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(
      coarsegraining_functions<scalar_type, k_dmn_t, q_dmn_t>& coarsegraining_functions_ref);

protected:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename K_dmn>
coarsegraining_routines<parameters_type, K_dmn>::coarsegraining_routines(parameters_type& parameters_ref)
    : parameters(parameters_ref), concurrency(parameters.get_concurrency()) {}

template <typename parameters_type, typename K_dmn>
void coarsegraining_routines<parameters_type, K_dmn>::compute_tetrahedron_mesh(int k_mesh_refinement,
                                                                               int number_of_periods) {
  math_algorithms::tetrahedron_mesh<typename K_dmn::parameter_type> mesh(k_mesh_refinement);

  quadrature_dmn::translate_according_to_period(number_of_periods, mesh);

  {
    coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_size() = 0;
    coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_weights().resize(0);
    coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_elements().resize(0);

    for (int l = 0; l < mesh.size(); l++)
      mesh[l].update_tetrahedron_domain(coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_size(),
                                        coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_weights(),
                                        coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_elements());

    coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_size() = 0;
    coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_weights().resize(0);
    coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements().resize(0);

    for (int l = 0; l < mesh.size(); l++)
      mesh[l].update_tetrahedron_domain(
          coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_size(),
          coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_weights(),
          coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements());
  }
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_routines<parameters_type, K_dmn>::compute_gaussian_mesh(
    int k_mesh_refinement, int gaussian_quadrature_rule, int number_of_periods) {
  {
    quadrature_dmn::initialize_Brillouin_zone(k_mesh_refinement, gaussian_quadrature_rule,
                                              number_of_periods);
  }

  {
    coarsegraining_domain<K_dmn, ORIGIN>::get_size() = quadrature_dmn::get_size();
    coarsegraining_domain<K_dmn, ORIGIN>::get_weights() = quadrature_dmn::get_weights();
    coarsegraining_domain<K_dmn, ORIGIN>::get_elements() = quadrature_dmn::get_elements();
  }

  {
    coarsegraining_domain<K_dmn, K>::get_size() = quadrature_dmn::get_size();
    coarsegraining_domain<K_dmn, K>::get_weights() = quadrature_dmn::get_weights();
    coarsegraining_domain<K_dmn, K>::get_elements() = quadrature_dmn::get_elements();
  }

  {
    coarsegraining_domain<K_dmn, K_PLUS_Q>::get_size() = quadrature_dmn::get_size();
    coarsegraining_domain<K_dmn, K_PLUS_Q>::get_weights() = quadrature_dmn::get_weights();
    coarsegraining_domain<K_dmn, K_PLUS_Q>::get_elements() = quadrature_dmn::get_elements();
  }

  {
    coarsegraining_domain<K_dmn, Q_MINUS_K>::get_size() = quadrature_dmn::get_size();
    coarsegraining_domain<K_dmn, Q_MINUS_K>::get_weights() = quadrature_dmn::get_weights();
    coarsegraining_domain<K_dmn, Q_MINUS_K>::get_elements() = quadrature_dmn::get_elements();
  }
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::wannier_interpolation(
    int K_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& f_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& f_q) {
  typedef interpolation_matrices<scalar_type, k_dmn_t, q_dmn_t> interpolation_matrices_type;

  if (not interpolation_matrices_type::is_initialized())
    interpolation_matrices_type::initialize(concurrency);

  LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& T = interpolation_matrices_type::get(K_ind);

  scalar_type alpha(1.);
  scalar_type beta(0.);

  scalar_type* A_ptr = &(reinterpret_cast<scalar_type(&)[2]>(f_k(0))[0]);  // &real(f_k(0));
  scalar_type* B_ptr = &T(0, 0);
  scalar_type* C_ptr = &(reinterpret_cast<scalar_type(&)[2]>(f_q(0))[0]);  // &real(f_q(0));

  int M = 2 * nu_nu::dmn_size();
  int K = k_dmn_t::dmn_size();
  int N = q_dmn_t::dmn_size();

  int LDA = 2 * nu_nu::dmn_size();
  int LDB = T.get_global_size().first;
  int LDC = 2 * nu_nu::dmn_size();

  dca::linalg::blas::gemm("N", "T", M, N, K, alpha, A_ptr, LDA, B_ptr, LDB, beta, C_ptr, LDC);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_I_q(
    tmp_scalar_type value,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q) {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        I_q(i, j, q_ind) = i == j ? value : 0.;
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_H_q(
    int K_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q) {
  wannier_interpolation(K_ind, H_k, H_q);
}

/*****************************************
 ***                                   ***
 ***         Routines for DCA          ***
 ***                                   ***
 *****************************************/

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q(
    int K_ind, int w_ind,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w>>& S_K_w,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q) {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        S_q(i, j, q_ind) = S_K_w(i, j, K_ind, w_ind);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    int K_ind, int w_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w>>& S_K,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(parameters.get_chemical_potential());
    i_wm_min_mu.imag(w::get_elements()[w_ind]);

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q(K_ind, w_ind, S_K, S_q);

  int nr_threads = parameters.get_nr_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_st(I_q, H_q, S_q,
                                                                                      G_q);
  else
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(
        nr_threads, I_q, H_q, S_q, G_q);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q(
    int K_ind, int w_ind,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w_REAL>>& S_K_w,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q) {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        S_q(i, j, q_ind) = S_K_w(i, j, K_ind, w_ind);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    int K_ind, int w_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_DCA, w_REAL>>& S_K,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(w_REAL::get_elements()[w_ind] + parameters.get_chemical_potential());
    i_wm_min_mu.imag(parameters.get_real_frequencies_off_set());

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q(K_ind, w_ind, S_K, S_q);

  int nr_threads = parameters.get_nr_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_st(I_q, H_q, S_q,
                                                                                      G_q);
  else
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(
        nr_threads, I_q, H_q, S_q, G_q);
}

/*!
 *
 * p 122 AGD
 *
 */
template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_t(
    int K_ind, int t_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  scalar_type f_val = 1;
  scalar_type t_val = t::get_elements()[t_ind];
  scalar_type beta = parameters.get_beta();

  f_val = t_val < 0 ? 1 : -1;
  t_val = t_val < 0 ? t_val + beta : t_val;

  compute_I_q(parameters.get_chemical_potential(), I_q);

  compute_H_q(K_ind, H_k, H_q);

  G_q = 0.;

  int nr_threads = parameters.get_nr_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_st(
        beta, f_val, t_val, I_q, H_q, G_q);
  else
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_t_mt(
        nr_threads, beta, f_val, t_val, I_q, H_q, G_q);
}

/*****************************************
 ***                                   ***
 ***         Routines for DCA+         ***
 ***                                   ***
 *****************************************/

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q_from_A_k(
    int K_ind, int w_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& A_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& A_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q) {
  wannier_interpolation(K_ind, A_k, A_q);

  scalar_type alpha = w::get_elements()[w_ind] > 0 ? 1 : -1;

  transform_to_alpha::backward(alpha, S_q, A_q);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    int K_ind, int w_ind, FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t>>& A_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& A_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t>>& G_q) {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(parameters.get_chemical_potential());
    i_wm_min_mu.imag(w::get_elements()[w_ind]);

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q_from_A_k(K_ind, w_ind, A_k, A_q, S_q);

  int nr_threads = parameters.get_nr_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_st(I_q, H_q, S_q,
                                                                                      G_q);
  else
    quadrature_integration<parameters_type, q_dmn_t>::quadrature_integration_G_q_w_mt(
        nr_threads, I_q, H_q, S_q, G_q);
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_ROUTINES_H
