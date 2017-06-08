// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides methods to compute the free Green's function on various domains.

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_FREE_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_FREE_GREENS_FUNCTION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/quadrature_integration.hpp"

namespace dca {
namespace phys {
// dca::phys::

// Computes the free Green's function G_0(\vec{k}, i\omega_n) for a dispersion relation
// \epsilon_\vec{k} that is a matrix in orbital-spin space.
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename MatsubaraFreqDmn>
void compute_G0_k_w(
    /*const*/ func::function<std::complex<Scalar>,
                             func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& eps_k,
    const Scalar mu,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>& G0_k_w) {
  const std::complex<Scalar> i(0, 1);  // complex i

  // Need to pass a zero self-energy function to quadrature_integration_G_q_w_st.
  /*const*/ func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>> zero;

  // Diagonal i \omega_n + \mu function.
  func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>
      i_omega_n_plus_mu;

  // Helper function to store the result for fixed frequency.
  func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>> g;

  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w) {
    // Compute diagonal i \omega_n + \mu function.
    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        i_omega_n_plus_mu(m, m, k) = i * MatsubaraFreqDmn::get_elements()[w] + mu;
      }
    }

    g = 0.;

    clustermapping::quadrature_integration<KDmn, OrbitalSpinDmn>::quadrature_integration_G_q_w_st(
        i_omega_n_plus_mu, eps_k, zero, g);

    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
        for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
          G0_k_w(m, n, k, w) = g(m, n, k);
        }
      }
    }
  }
}

// Computes the free Green's function G_0(\vec{k}, \tau) for a dispersion relation \epsilon_\vec{k}
// that is a matrix in orbital-spin space.
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename ImagTimeDmn>
void compute_G0_k_t(
    /*const*/ func::function<std::complex<Scalar>,
                             func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& eps_k,
    const Scalar mu, const Scalar beta,
    func::function<Scalar, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, ImagTimeDmn>>& G0_k_t) {
  // Diagonal \mu function.
  func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>> mu_function;
  for (int k = 0; k < KDmn::dmn_size(); ++k) {
    for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
      mu_function(m, m, k) = mu;
    }
  }

  // Helper function to store the result for fixed time.
  func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>> g;

  for (int t = 0; t < ImagTimeDmn::dmn_size(); ++t) {
    Scalar tau = ImagTimeDmn::get_elements()[t];
    Scalar sign = -1;  // quadrature_integration_G_q_t_st requires Scalar.

    // G_0(\tau) = -G_0(\tau+\beta)
    if (tau < 0) {
      tau += beta;
      sign = 1;
    }

    g = 0.;

    clustermapping::quadrature_integration<KDmn, OrbitalSpinDmn>::quadrature_integration_G_q_t_st(
        beta, sign, tau, mu_function, eps_k, g);

    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
        for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
          if (g(m, n, k).imag() > 1.e-6)
            throw std::logic_error("G_0(\vec{k}, \tau) is real!");
          G0_k_t(m, n, k, t) = g(m, n, k).real();
        }
      }
    }
  }
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_FREE_GREENS_FUNCTION_HPP
