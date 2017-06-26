// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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

namespace dca {
namespace phys {
// dca::phys::

// Computes the interacting Green's function G(\vec{k}, i\omega_n) from the non-interacting Green's
// function G_0(\vec{k}, i\omega_n) = [i\omega_n + \mu - H_0(\vec{k})]^{-1} and the self-energy
// \Sigma(\vec{k}, i\omega_n) via the Dyson equation,
// G = [G_0^{-1} - \Sigma]^{-1}.
template <typename Scalar, typename OrbitalSpinDmn, typename KDmn, typename MatsubaraFreqDmn,
          typename ConcurrencyType>
void compute_G_k_w(
    /*const*/ func::function<std::complex<Scalar>,
                             func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn>>& H0_k,
    /*const*/ func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn,
                                                                      KDmn, MatsubaraFreqDmn>>& S_k_w,
    const Scalar mu, const ConcurrencyType& concurrency,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>& G_k_w) {
  // Work space for inverse.
  linalg::Matrix<std::complex<Scalar>, linalg::CPU> G_inv("G_inv", OrbitalSpinDmn::dmn_size());
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<std::complex<Scalar>, linalg::CPU> work;

  const std::complex<Scalar> i(0., 1.);

  // Distribute the work amongst the processes.
  func::dmn_variadic<KDmn, MatsubaraFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency.get_bounds(k_w_dmn_obj);
  int coor[2];

  G_k_w = 0.;

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, coor);
    const auto k_ind = coor[0];
    const auto w_ind = coor[1];
    const auto w_val = MatsubaraFreqDmn::get_elements()[w_ind];

    // Compute G^{-1} for fixed k-vector and Matsubara frequency.
    for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        G_inv(m, n) = -H0_k(m, n, k_ind) - S_k_w(m, n, k_ind, w_ind);
        if (m == n)
          G_inv(m, n) += i * w_val + mu;
      }
    }

    linalg::matrixop::inverse(G_inv, ipiv, work);

    for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n)
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m)
        G_k_w(m, n, k_ind, w_ind) = G_inv(m, n);
  }

  concurrency.sum(G_k_w);
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_GREENS_FUNCTION_HPP
