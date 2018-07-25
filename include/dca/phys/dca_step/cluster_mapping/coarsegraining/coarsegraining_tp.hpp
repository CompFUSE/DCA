// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class performs the coarsegraining of two-particle functions.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_TP_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_TP_HPP

#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/geometry/gaussian_quadrature/gaussian_quadrature_domain.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_mesh.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_routines.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/interpolation_matrices.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/transform_to_alpha.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/util/plot.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type, typename K_dmn>
class coarsegraining_tp : public coarsegraining_routines<parameters_type> {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using k_cluster_type = typename K_dmn::parameter_type;

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;
  using host_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_SP,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_HOST = func::dmn_0<host_k_cluster_type>;

  using scalar_type = double;
  using complex_type = std::complex<scalar_type>;

  using tetrahedron_dmn = func::dmn_0<math::geometry::tetrahedron_mesh<K_dmn>>;
  using quadrature_dmn = math::geometry::gaussian_quadrature_domain<tetrahedron_dmn>;

  using q_dmn = func::dmn_0<coarsegraining_domain<K_dmn, K>>;
  using q_plus_Q_dmn = func::dmn_0<coarsegraining_domain<K_dmn, K_PLUS_Q>>;
  using Q_min_q_dmn = func::dmn_0<coarsegraining_domain<K_dmn, Q_MINUS_K>>;

  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using b_b = func::dmn_variadic<b, b>;

  using nu_nu_q = func::dmn_variadic<nu, nu, q_dmn>;
  using nu_nu_q_plus_Q = func::dmn_variadic<nu, nu, q_plus_Q_dmn>;
  using nu_nu_Q_min_q = func::dmn_variadic<nu, nu, Q_min_q_dmn>;

  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

public:
  coarsegraining_tp(parameters_type& parameters_ref);

  // DCA coarsegraining
  template <typename w_dmn_t>
  void execute(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& Sigma,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi);

  // DCA+ coarsegraining
  template <typename w_dmn_t>
  void execute(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& Sigma,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi);

  template <typename w_dmn_t>
  void plot(func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi);

private:
  // DCA coarsegraining
  template <typename w_dmn_t>
  void compute_tp(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& Sigma,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi);

  // DCA+ coarsegraining
  template <typename w_dmn_t>
  void compute_tp(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& Sigma,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi);

  // DCA coarsegraining
  template <typename w_dmn_t>
  void compute_phi(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& phi);

  // DCA+ coarsegraining
  template <typename w_dmn_t>
  void compute_phi(
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& S_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& phi);

  void find_w1_and_w2(std::vector<double>& elements, int& w_ind, int& w1, int& w2);

  void compute_bubble(
      func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, q_dmn>>& bubble);

  double get_integration_factor();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  func::function<scalar_type, q_dmn> w_q;

  func::function<std::complex<scalar_type>, nu_nu_q> I_q;
  func::function<std::complex<scalar_type>, nu_nu_q> H_q;
  func::function<std::complex<scalar_type>, nu_nu_q> S_q;
  func::function<std::complex<scalar_type>, nu_nu_q> A_q;
  func::function<std::complex<scalar_type>, nu_nu_q> G_q;

  func::function<std::complex<scalar_type>, nu_nu_q_plus_Q> I_q_plus_Q;
  func::function<std::complex<scalar_type>, nu_nu_q_plus_Q> H_q_plus_Q;
  func::function<std::complex<scalar_type>, nu_nu_q_plus_Q> S_q_plus_Q;
  func::function<std::complex<scalar_type>, nu_nu_q_plus_Q> A_q_plus_Q;
  func::function<std::complex<scalar_type>, nu_nu_q_plus_Q> G_q_plus_Q;

  func::function<std::complex<scalar_type>, nu_nu_Q_min_q> I_Q_min_q;
  func::function<std::complex<scalar_type>, nu_nu_Q_min_q> H_Q_min_q;
  func::function<std::complex<scalar_type>, nu_nu_Q_min_q> S_Q_min_q;
  func::function<std::complex<scalar_type>, nu_nu_Q_min_q> A_Q_min_q;
  func::function<std::complex<scalar_type>, nu_nu_Q_min_q> G_Q_min_q;

  func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, q_dmn>> bubble_q;
};

template <typename parameters_type, typename K_dmn>
coarsegraining_tp<parameters_type, K_dmn>::coarsegraining_tp(parameters_type& parameters_ref)
    : coarsegraining_routines<parameters_type>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      w_q("w_q"),

      I_q("I_q"),
      H_q("H_q"),
      S_q("S_q"),
      A_q("A_q"),
      G_q("G_q"),

      I_q_plus_Q("I_q_plus_Q"),
      H_q_plus_Q("H_q_plus_Q"),
      S_q_plus_Q("S_q_plus_Q"),
      A_q_plus_Q("A_q_plus_Q"),
      G_q_plus_Q("G_q_plus_Q"),

      I_Q_min_q("I_Q_min_q"),
      H_Q_min_q("H_Q_min_q"),
      S_Q_min_q("S_Q_min_q"),
      A_Q_min_q("A_Q_min_q"),
      G_Q_min_q("G_Q_min_q"),

      bubble_q("bubble_q") {
  for (int l = 0; l < w_q.size(); l++)
    w_q(l) = quadrature_dmn::get_weights()[l];
}

// DCA-coarsegraining
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::execute(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& Sigma,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi) {
  int Q_ind = domains::cluster_operations::index(parameters.get_four_point_momentum_transfer(),
                                                 K_dmn::get_elements(), K_dmn::parameter_type::SHAPE);

  switch (parameters.get_four_point_type()) {
    case PARTICLE_HOLE_CHARGE:
    case PARTICLE_HOLE_MAGNETIC:
    case PARTICLE_HOLE_TRANSVERSE: {
      interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);
      interpolation_matrices<scalar_type, k_HOST, q_plus_Q_dmn>::initialize(concurrency, Q_ind);

      compute_tp(H_k, Sigma, chi);
    } break;

    case PARTICLE_PARTICLE_UP_DOWN: {
      interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);
      interpolation_matrices<scalar_type, k_HOST, Q_min_q_dmn>::initialize(concurrency, Q_ind);

      compute_phi(H_k, Sigma, chi);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

// DCA+-coarsegraining
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::execute(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& Sigma,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi) {
  int Q_ind = domains::cluster_operations::index(parameters.get_four_point_momentum_transfer(),
                                                 K_dmn::get_elements(), K_dmn::parameter_type::SHAPE);

  switch (parameters.get_four_point_type()) {
    case PARTICLE_HOLE_CHARGE:
    case PARTICLE_HOLE_MAGNETIC:
    case PARTICLE_HOLE_TRANSVERSE: {
      interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);
      interpolation_matrices<scalar_type, k_HOST, q_plus_Q_dmn>::initialize(concurrency, Q_ind);

      compute_tp(H_k, Sigma, chi);
    } break;

    case PARTICLE_PARTICLE_UP_DOWN: {
      interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);
      interpolation_matrices<scalar_type, k_HOST, Q_min_q_dmn>::initialize(concurrency, Q_ind);

      compute_phi(H_k, Sigma, chi);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::plot(
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi) {
  {
    util::Plot plot;

    func::function<scalar_type, w_dmn_t> phi_w("phi_w");

    for (int m2 = 0; m2 < b::dmn_size(); m2++) {
      for (int m1 = 0; m1 < b::dmn_size(); m1++) {
        for (int n2 = 0; n2 < b::dmn_size(); n2++) {
          for (int n1 = 0; n1 < b::dmn_size(); n1++) {
            for (int k_ind = 0; k_ind < K_dmn::dmn_size(); k_ind++) {
              if (abs(K_dmn::get_elements()[k_ind][0] + K_dmn::get_elements()[k_ind][1] - M_PI) <
                  1.e-6) {
                for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++)
                  phi_w(w_ind) = real(chi(n1, n2, m1, m2, k_ind, w_ind));

                plot.plotLinesPoints(phi_w);
              }
            }
          }
        }
      }
    }
  }

  {
    util::Plot plot;

    func::function<scalar_type, w_dmn_t> phi_w("phi_w");

    for (int m2 = 0; m2 < b::dmn_size(); m2++) {
      for (int m1 = 0; m1 < b::dmn_size(); m1++) {
        for (int n2 = 0; n2 < b::dmn_size(); n2++) {
          for (int n1 = 0; n1 < b::dmn_size(); n1++) {
            for (int k_ind = 0; k_ind < K_dmn::dmn_size(); k_ind++) {
              if (abs(K_dmn::get_elements()[k_ind][0] + K_dmn::get_elements()[k_ind][1] - M_PI) <
                  1.e-6) {
                for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++)
                  phi_w(w_ind) = imag(chi(n1, n2, m1, m2, k_ind, w_ind));

                plot.plot(phi_w);
              }
            }
          }
        }
      }
    }
  }

  {
    std::vector<double> x(0);
    std::vector<double> y(0);
    std::vector<double> z(0);

    for (int k_ind = 0; k_ind < K_dmn::dmn_size(); k_ind++) {
      x.push_back(K_dmn::get_elements()[k_ind][0]);
      y.push_back(K_dmn::get_elements()[k_ind][1]);
      z.push_back(real(phi(0, 0, 0, 0, k_ind, w_dmn_t::dmn_size() / 2)));
    }

    util::Plot::heatMap(x, y, z);
  }
}

// DCA coarsegaining where K-dmn is the cluster-domain
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::compute_tp(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi) {
  assert(k_DCA::get_elements() == K_dmn::get_elements());

  chi = 0.;

  K_dmn k_domain;
  std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

  // S_K_plus_Q_w(K) = S_K_w(K+Q)
  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>> S_K_plus_Q_w;

  int Q_ind = domains::cluster_operations::index(parameters.get_four_point_momentum_transfer(),
                                                 k_DCA::get_elements(), k_DCA::parameter_type::SHAPE);

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
      int K_plus_Q_ind = k_DCA::parameter_type::add(k_ind, Q_ind);

      for (int nu_2 = 0; nu_2 < nu::dmn_size(); ++nu_2) {
        for (int nu_1 = 0; nu_1 < nu::dmn_size(); ++nu_1) {
          S_K_plus_Q_w(nu_1, nu_2, k_ind, w_ind) = S_K_w(nu_1, nu_2, K_plus_Q_ind, w_ind);
        }
      }
    }
  }

  for (int k_ind = bounds.first; k_ind < bounds.second; k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      int w_1, w_2;

      find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

      {
        coarsegraining_routines<parameters_type>::compute_G_q_w(k_ind, w_1, H_k, S_K_w, I_q,
                                                                       H_q, S_q, G_q);
      }

      {
        coarsegraining_routines<parameters_type>::compute_G_q_w(
            k_ind, w_2, H_k, S_K_plus_Q_w, I_q_plus_Q, H_q_plus_Q, S_q_plus_Q, G_q_plus_Q);
      }

      compute_bubble(bubble_q);

      {
        double factor = get_integration_factor();

        for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
          for (int n1 = 0; n1 < b::dmn_size(); n1++)
            for (int n2 = 0; n2 < b::dmn_size(); n2++)
              for (int m1 = 0; m1 < b::dmn_size(); m1++)
                for (int m2 = 0; m2 < b::dmn_size(); m2++)
                  chi(n1, n2, m1, m2, k_ind, w_ind) +=
                      factor * w_q(q_ind) * bubble_q(n1, n2, m1, m2, q_ind);
      }
    }
  }

  concurrency.sum(chi);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    chi /= V_K;
  }
}

// DCA+ coarsegaining where K-dmn is the host
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::compute_tp(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& S_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& chi) {
  chi = 0.;

  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>> A_k("A_k");
  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>> A_k_w("A_k_w");

  latticemapping::transform_to_alpha::forward(1., S_k_w, A_k_w);

  K_dmn k_domain;
  std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

  for (int k_ind = bounds.first; k_ind < bounds.second; k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      int w_1, w_2;

      find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

      {
        for (int k = 0; k < k_HOST::dmn_size(); k++)
          for (int j = 0; j < nu::dmn_size(); j++)
            for (int i = 0; i < nu::dmn_size(); i++)
              A_k(i, j, k) = A_k_w(i, j, k, w_1);

        coarsegraining_routines<parameters_type>::compute_G_q_w(k_ind, w_1, H_k, A_k, I_q,
                                                                       H_q, A_q, S_q, G_q);
      }

      {
        for (int k = 0; k < k_HOST::dmn_size(); k++)
          for (int j = 0; j < nu::dmn_size(); j++)
            for (int i = 0; i < nu::dmn_size(); i++)
              A_k(i, j, k) = A_k_w(i, j, k, w_2);

        coarsegraining_routines<parameters_type>::compute_G_q_w(
            k_ind, w_2, H_k, A_k, I_q_plus_Q, H_q_plus_Q, A_q_plus_Q, S_q_plus_Q, G_q_plus_Q);
      }

      compute_bubble(bubble_q);

      {
        double factor = get_integration_factor();

        for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
          for (int n1 = 0; n1 < b::dmn_size(); n1++)
            for (int n2 = 0; n2 < b::dmn_size(); n2++)
              for (int m1 = 0; m1 < b::dmn_size(); m1++)
                for (int m2 = 0; m2 < b::dmn_size(); m2++)
                  chi(n1, n2, m1, m2, k_ind, w_ind) +=
                      factor * w_q(q_ind) * bubble_q(n1, n2, m1, m2, q_ind);
      }
    }
  }

  concurrency.sum(chi);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    chi /= V_K;
  }
}

// DCA coarsegaining where K-dmn is the cluster-domain
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::compute_phi(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>>& S_K_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& phi) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t start " << __FUNCTION__ << " ... " << dca::util::print_time();

  assert(k_DCA::get_elements() == K_dmn::get_elements());

  phi = 0.;

  K_dmn k_domain;
  std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

  // S_Q_min_K_w(K) = S_K_w(Q-K)
  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_DCA, w>> S_Q_min_K_w;

  int Q_ind = domains::cluster_operations::index(parameters.get_four_point_momentum_transfer(),
                                                 k_DCA::get_elements(), k_DCA::parameter_type::SHAPE);

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
      int Q_min_K_ind = k_DCA::parameter_type::subtract(k_ind, Q_ind);

      for (int nu_2 = 0; nu_2 < nu::dmn_size(); ++nu_2) {
        for (int nu_1 = 0; nu_1 < nu::dmn_size(); ++nu_1) {
          S_Q_min_K_w(nu_1, nu_2, k_ind, w_ind) = S_K_w(nu_1, nu_2, Q_min_K_ind, w_ind);
        }
      }
    }
  }

  for (int k_ind = bounds.first; k_ind < bounds.second; k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      int w_1, w_2;

      find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

      {
        coarsegraining_routines<parameters_type>::compute_G_q_w(k_ind, w_1, H_k, S_K_w, I_q,
                                                                       H_q, S_q, G_q);
      }

      {
        coarsegraining_routines<parameters_type>::compute_G_q_w(
            k_ind, w_2, H_k, S_Q_min_K_w, I_Q_min_q, H_Q_min_q, S_Q_min_q, G_Q_min_q);
      }

      compute_bubble(bubble_q);

      {
        double factor = get_integration_factor();

        for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
          for (int m2 = 0; m2 < b::dmn_size(); m2++)
            for (int m1 = 0; m1 < b::dmn_size(); m1++)
              for (int n2 = 0; n2 < b::dmn_size(); n2++)
                for (int n1 = 0; n1 < b::dmn_size(); n1++)
                  phi(n1, n2, m1, m2, k_ind, w_ind) +=
                      factor * w_q(q_ind) * bubble_q(n1, n2, m1, m2, q_ind);
      }
    }
  }

  concurrency.sum(phi);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    phi /= V_K;
  }
}

// DCA+ coarsegaining where K-dmn is the host
template <typename parameters_type, typename K_dmn>
template <typename w_dmn_t>
void coarsegraining_tp<parameters_type, K_dmn>::compute_phi(
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>>& H_k,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>>& S_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, K_dmn, w_dmn_t>>& phi) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t start " << __FUNCTION__ << " ... " << dca::util::print_time();

  phi = 0.;

  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST>> A_k("A_k");
  func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_HOST, w>> A_k_w("A_k_w");

  latticemapping::transform_to_alpha::forward(1., S_k_w, A_k_w);

  K_dmn k_domain;
  std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

  for (int k_ind = bounds.first; k_ind < bounds.second; k_ind++) {
    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
      int w_1, w_2;

      find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

      {
        for (int k = 0; k < k_HOST::dmn_size(); k++)
          for (int j = 0; j < nu::dmn_size(); j++)
            for (int i = 0; i < nu::dmn_size(); i++)
              A_k(i, j, k) = A_k_w(i, j, k, w_1);

        coarsegraining_routines<parameters_type>::compute_G_q_w(k_ind, w_1, H_k, A_k, I_q,
                                                                       H_q, A_q, S_q, G_q);
      }

      {
        for (int k = 0; k < k_HOST::dmn_size(); k++)
          for (int j = 0; j < nu::dmn_size(); j++)
            for (int i = 0; i < nu::dmn_size(); i++)
              A_k(i, j, k) = A_k_w(i, j, k, w_2);

        coarsegraining_routines<parameters_type>::compute_G_q_w(
            k_ind, w_2, H_k, A_k, I_Q_min_q, H_Q_min_q, A_Q_min_q, S_Q_min_q, G_Q_min_q);
      }

      compute_bubble(bubble_q);

      {
        double factor = get_integration_factor();

        for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
          for (int m2 = 0; m2 < b::dmn_size(); m2++)
            for (int m1 = 0; m1 < b::dmn_size(); m1++)
              for (int n2 = 0; n2 < b::dmn_size(); n2++)
                for (int n1 = 0; n1 < b::dmn_size(); n1++)
                  phi(n1, n2, m1, m2, k_ind, w_ind) +=
                      factor * w_q(q_ind) * bubble_q(n1, n2, m1, m2, q_ind);
      }
    }
  }

  concurrency.sum(phi);

  {
    scalar_type V_K = 0;
    for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++)
      V_K += w_q(q_ind);

    phi /= V_K;
  }

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t end  ... " << dca::util::print_time();
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_tp<parameters_type, K_dmn>::find_w1_and_w2(std::vector<double>& elements,
                                                               int& w_ind, int& w1, int& w2) {
  int W_ind = parameters.get_four_point_frequency_transfer();

  for (int l = 0; l < w::dmn_size(); l++)
    if (std::abs(elements[w_ind] - w::get_elements()[l]) < 1.e-6)
      w1 = l;

  assert(std::abs(w::get_elements()[w1] - elements[w_ind]) < 1.e-6);

  switch (parameters.get_four_point_type()) {
    case PARTICLE_HOLE_CHARGE:
    case PARTICLE_HOLE_MAGNETIC:
    case PARTICLE_HOLE_TRANSVERSE: {
      w2 = w1 + W_ind;
    } break;

    case PARTICLE_PARTICLE_UP_DOWN: {
      w2 = W_ind + (w::dmn_size() - 1 - w1);
      assert(std::abs(w::get_elements()[w1] + w::get_elements()[w::dmn_size() - 1 - w1]) < 1.e-6);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename parameters_type, typename K_dmn>
void coarsegraining_tp<parameters_type, K_dmn>::compute_bubble(
    func::function<std::complex<scalar_type>, func::dmn_variadic<b_b, b_b, q_dmn>>& bubble) {
  bubble = 0.;

  for (int q_ind = 0; q_ind < q_dmn::dmn_size(); q_ind++) {
    for (int n1 = 0; n1 < b::dmn_size(); n1++) {
      for (int n2 = 0; n2 < b::dmn_size(); n2++) {
        for (int m1 = 0; m1 < b::dmn_size(); m1++) {
          for (int m2 = 0; m2 < b::dmn_size(); m2++) {
            switch (parameters.get_four_point_type()) {
              case PARTICLE_HOLE_TRANSVERSE:
                bubble(n1, n2, m1, m2, q_ind) +=
                    G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                break;

              case PARTICLE_HOLE_MAGNETIC:
                bubble(n1, n2, m1, m2, q_ind) +=
                    G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                break;

              case PARTICLE_HOLE_CHARGE:
                bubble(n1, n2, m1, m2, q_ind) +=
                    G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                break;

              case PARTICLE_PARTICLE_UP_DOWN:
                bubble(n1, n2, m1, m2, q_ind) +=
                    G_q(n1, e_UP, m1, e_UP, q_ind) * G_Q_min_q(n2, e_UP, m2, e_UP, q_ind);
                break;

              default:
                throw std::logic_error(__FUNCTION__);
            }
          }
        }
      }
    }
  }
}

template <typename parameters_type, typename K_dmn>
double coarsegraining_tp<parameters_type, K_dmn>::get_integration_factor() {
  switch (parameters.get_four_point_type()) {
    case PARTICLE_HOLE_TRANSVERSE:
      return -1.;
      break;

    case PARTICLE_HOLE_MAGNETIC:
      return -1.;
      break;

    case PARTICLE_HOLE_CHARGE:
      return -2.;
      break;

    case PARTICLE_PARTICLE_UP_DOWN:
      return 1;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_TP_HPP
