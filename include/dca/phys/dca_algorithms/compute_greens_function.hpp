// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a method to compute the (interacting) Green's function G(\vec{k}, i\omega_n).

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_GREENS_FUNCTION_HPP

#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/parallel/util/get_bounds.hpp"

namespace dca {
namespace phys {
// dca::phys::

// Computes the interacting Green's function G(\vec{k}, i\omega_n) from the non-interacting Green's
// function G_0(\vec{k}, i\omega_n) = [i\omega_n + \mu - H_0(\vec{k})]^{-1} and the self-energy
// \Sigma(\vec{k}, i\omega_n) via the Dyson equation,
// G = [G_0^{-1} - \Sigma]^{-1}.
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename MatsubaraFreqDmn>
void compute_G_k_w(
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& H0_k,
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>& S_k_w,
    const Scalar mu, const int n_threads,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>& G_k_w) {
  G_k_w = 0.;

  parallel::stdthread().execute(n_threads, [&](int id, int n_threads) {
    // Work space for inverse.
    linalg::Matrix<std::complex<Scalar>, linalg::CPU> G_inv("G_inv", OrbitalSpinDmn::dmn_size());
    linalg::Vector<int, linalg::CPU> ipiv;
    linalg::Vector<std::complex<Scalar>, linalg::CPU> work;

    const std::complex<Scalar> i(0., 1.);

    // Distribute the work amongst the threads.
    const auto bounds = parallel::util::getBounds(id, n_threads, MatsubaraFreqDmn());

    for (int w_ind = bounds.first; w_ind < bounds.second; ++w_ind)
      for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
        const auto w_val = MatsubaraFreqDmn::get_elements()[w_ind];

        // Compute G^{-1} for fixed k-vector and Matsubara frequency.
        for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
          for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
            G_inv(m, n) = -H0_k(m, n, k_ind) - S_k_w(m, n, k_ind, w_ind);
            if (m == n)
              G_inv(m, n) += i * w_val + mu;
          }
        }

        linalg::matrixop::smallInverse(G_inv, ipiv, work);

        for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n)
          for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m)
            G_k_w(m, n, k_ind, w_ind) = G_inv(m, n);
      }
  });
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_GREENS_FUNCTION_HPP
