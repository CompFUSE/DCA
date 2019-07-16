// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class provides routines for the coarsegraining.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_ROUTINES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_ROUTINES_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/geometry/gaussian_quadrature/gaussian_quadrature_domain.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_mesh.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_domain.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/interpolation_matrices.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/quadrature_integration.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/transform_to_alpha.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type,
          typename K_dmn = typename ClusterDomainAliases<parameters_type::lattice_dimension>::KClusterDmn>
class coarsegraining_routines {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using k_cluster_type = typename K_dmn::parameter_type;

  using tetrahedron_dmn = func::dmn_0<math::geometry::tetrahedron_mesh<K_dmn>>;
  using quadrature_dmn = math::geometry::gaussian_quadrature_domain<tetrahedron_dmn>;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;

  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

public:
  coarsegraining_routines(parameters_type& parameters_ref);

  // Transforms f_k into f_q by interpolating the shifted inverse as
  // f_q = [interpolation((f_q - i*alpha)^-1)]^-1 + i*alpha.
  // In:  f_k.
  // Out: f_q.
  // Postcondition: f_k is shifted and inverted.
  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void wannierInterpolationWithAlphaTransform(
      int K_ind, double alpha,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& f_q) const;

protected:
  void compute_tetrahedron_mesh(int k_mesh_refinement, int number_of_periods) const;

  void compute_gaussian_mesh(int k_mesh_refinement, int gaussian_quadrature_rule,
                             int number_of_periods) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void wannier_interpolation(
      int K_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& f_q) const;

  template <typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
  void compute_I_q(
      tmp_scalar_type value,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_H_q(
      int K_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q) const;

  /***********************************************************************************
   ***
   ***         Routines for DCA
   ***
   ***********************************************************************************/

  template <typename scalar_type, typename q_dmn_t>
  void compute_S_q(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const;

  template <typename scalar_type, typename q_dmn_t>
  void compute_S_q(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>>& S_K_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>>& S_K,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_t(
      int K_ind, int t_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const;

  /***********************************************************************************
   ***                                                                             ***
   ***         Routines for DCA+
   ***                                                                             ***
   ***********************************************************************************/

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_S_q_from_A_k(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& A_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& A_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const;

  template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void compute_G_q_w(
      int K_ind, int w_ind,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
      const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& A_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& A_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const;

protected:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename K_dmn>
coarsegraining_routines<parameters_type, K_dmn>::coarsegraining_routines(parameters_type& parameters_ref)
    : parameters(parameters_ref), concurrency(parameters.get_concurrency()) {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    linalg::lapack::silenceLapack();

    compute_tetrahedron_mesh(parameters.get_k_mesh_recursion(),
                             parameters.get_coarsegraining_periods());

    compute_gaussian_mesh(parameters.get_k_mesh_recursion(), parameters.get_quadrature_rule(),
                          parameters.get_coarsegraining_periods());
  });
}

template <class Parameters, typename K_dmn>
template <typename ScalarType, typename KDmn, typename QDmn>
void coarsegraining_routines<Parameters, K_dmn>::wannierInterpolationWithAlphaTransform(
    const int k_ind, const double alpha,
    func::function<std::complex<ScalarType>, func::dmn_variadic<nu, nu, KDmn>>& f_k,
    func::function<std::complex<ScalarType>, func::dmn_variadic<nu, nu, QDmn>>& f_q) const {
  latticemapping::transform_to_alpha::forward(alpha, f_k);
  wannier_interpolation(k_ind, f_k, f_q);
  latticemapping::transform_to_alpha::backward(alpha, f_q);
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_routines<parameters_type, K_dmn>::compute_tetrahedron_mesh(
    const int k_mesh_refinement, const int number_of_periods) const {
  math::geometry::tetrahedron_mesh<typename K_dmn::parameter_type> mesh(k_mesh_refinement);

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
    const int k_mesh_refinement, const int gaussian_quadrature_rule,
    const int number_of_periods) const {
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
    const int K_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& f_q) const {
  typedef interpolation_matrices<scalar_type, k_dmn_t, q_dmn_t> interpolation_matrices_type;

  if (!interpolation_matrices_type::is_initialized())
    interpolation_matrices_type::initialize(concurrency);

  const auto T = interpolation_matrices_type::get(K_ind);

  scalar_type alpha(1.);
  scalar_type beta(0.);

  // Interpolate real, imaginary part and spin bands independently.
  const scalar_type* A_ptr = reinterpret_cast<const scalar_type*>(f_k.values());
  const scalar_type* B_ptr = &T(0, 0);
  scalar_type* C_ptr = reinterpret_cast<scalar_type*>(f_q.values());

  int M = 2 * nu_nu::dmn_size();
  int K = k_dmn_t::dmn_size();
  int N = q_dmn_t::dmn_size();

  int LDA = 2 * nu_nu::dmn_size();
  int LDB = T.leadingDimension();
  int LDC = 2 * nu_nu::dmn_size();

  dca::linalg::blas::gemm("N", "T", M, N, K, alpha, A_ptr, LDA, B_ptr, LDB, beta, C_ptr, LDC);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_I_q(
    const tmp_scalar_type value,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q) const {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        I_q(i, j, q_ind) = i == j ? value : 0.;
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_H_q(
    const int K_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q) const {
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
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        S_q(i, j, q_ind) = S_K_w(i, j, K_ind, w_ind);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(parameters.get_chemical_potential());
    i_wm_min_mu.imag(w::get_elements()[w_ind]);

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q(K_ind, w_ind, S_K, S_q);

  int nr_threads = parameters.get_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_w_st(I_q, H_q, S_q, G_q);
  else
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_w_mt(nr_threads, I_q, H_q, S_q, G_q);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q(
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>>& S_K_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const {
  for (int q_ind = 0; q_ind < q_dmn_t::dmn_size(); q_ind++)
    for (int j = 0; j < nu::dmn_size(); j++)
      for (int i = 0; i < nu::dmn_size(); i++)
        S_q(i, j, q_ind) = S_K_w(i, j, K_ind, w_ind);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w_REAL>>& S_K,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(w_REAL::get_elements()[w_ind] + parameters.get_chemical_potential());
    i_wm_min_mu.imag(parameters.get_imaginary_damping());

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q(K_ind, w_ind, S_K, S_q);

  int nr_threads = parameters.get_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<q_dmn_t, nu, typename parameters_type::ThreadingType>::quadrature_integration_G_q_w_st(
        I_q, H_q, S_q, G_q);
  else
    quadrature_integration<q_dmn_t, nu, typename parameters_type::ThreadingType>::quadrature_integration_G_q_w_mt(
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
    const int K_ind, const int t_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const {
  scalar_type f_val = 1;
  scalar_type t_val = t::get_elements()[t_ind];
  scalar_type beta = parameters.get_beta();

  f_val = t_val < 0 ? 1 : -1;
  t_val = t_val < 0 ? t_val + beta : t_val;

  compute_I_q(parameters.get_chemical_potential(), I_q);

  compute_H_q(K_ind, H_k, H_q);

  G_q = 0.;

  int nr_threads = parameters.get_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_t_st(beta, f_val, t_val, I_q, H_q, G_q);
  else
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_t_mt(nr_threads, beta, f_val, t_val, I_q, H_q, G_q);
}

/*****************************************
 ***                                   ***
 ***         Routines for DCA+         ***
 ***                                   ***
 *****************************************/

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q_from_A_k(
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& A_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& A_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q) const {
  wannier_interpolation(K_ind, A_k, A_q);

  scalar_type alpha = w::get_elements()[w_ind] > 0 ? 1 : -1;

  latticemapping::transform_to_alpha::backward(alpha, S_q, A_q);
}

template <typename parameters_type, typename K_dmn>
template <typename scalar_type, typename k_dmn_t, typename q_dmn_t>
void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(
    const int K_ind, const int w_ind,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& H_k,
    const func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& A_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& I_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& H_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& A_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& S_q,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, q_dmn_t>>& G_q) const {
  {
    std::complex<scalar_type> i_wm_min_mu;

    i_wm_min_mu.real(parameters.get_chemical_potential());
    i_wm_min_mu.imag(w::get_elements()[w_ind]);

    compute_I_q(i_wm_min_mu, I_q);
  }

  compute_H_q(K_ind, H_k, H_q);

  compute_S_q_from_A_k(K_ind, w_ind, A_k, A_q, S_q);

  int nr_threads = parameters.get_coarsegraining_threads();

  if (nr_threads == 1)
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_w_st(I_q, H_q, S_q, G_q);
  else
    quadrature_integration<scalar_type, q_dmn_t, nu, typename parameters_type::ThreadingType>::
        quadrature_integration_G_q_w_mt(nr_threads, I_q, H_q, S_q, G_q);
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_ROUTINES_HPP
