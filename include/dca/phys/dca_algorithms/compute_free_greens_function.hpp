// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides methods to compute the free Green's function on various domains.

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_FREE_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_FREE_GREENS_FUNCTION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_algorithms/compute_greens_function.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/quadrature_integration.hpp"

namespace dca {
namespace phys {
// dca::phys::

// Computes the free Matsubara Green's function G_0(\vec{k}, i\omega_n) from the non-interacting
// Hamiltonian H_0(\vec{k}).
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename MatsubaraFreqDmn,
          typename ConcurrencyType>
void compute_G0_k_w(
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& H0_k,
    const Scalar mu, const ConcurrencyType& concurrency,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>& G0_k_w) {
  // Call compute_G_k_w with vanishing self-energy.
  const func::function<std::complex<Scalar>,
                       func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>
      zero;
  compute_G_k_w(H0_k, zero, mu, concurrency, G0_k_w);
}

// Computes the free imaginary time Green's function G_0(\vec{k}, \tau) from the non-interacting
// Hamiltonian H_0(\vec{k}).
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename ImagTimeDmn>
void compute_G0_k_t(
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& H0_k,
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
        beta, sign, tau, mu_function, H0_k, g);

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
