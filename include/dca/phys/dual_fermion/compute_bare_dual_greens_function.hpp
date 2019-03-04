// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a method to compute the bare dual Green's function,
//     \tilde{G}^{(0)}(\vec{K}, \vec{K'}, \tilde{\vec{k}}, \omega_n)
//         = - G^{(c)}(\vec{K}, \omega_n)
//           * [G^{(c)}(\vec{K'}, \omega_n)
//              + ((\Delta(\vec{K'}, \omega_n) + \bar{\vareps}_\vec{K}) \delta_{\vec{K}, \vec{K'}}
//                 - \hat{\vareps}}_{\vec{K}, \vec{K'}, \tilde{\vec{k}}})^{-1}]^{-1}
//           * G^{(c)}(\vec{K'}, \omega_n) .

#ifndef DCA_PHYS_DUAL_FERMION_COMPUTE_BARE_DUAL_GREENS_FUNCTION_HPP
#define DCA_PHYS_DUAL_FERMION_COMPUTE_BARE_DUAL_GREENS_FUNCTION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Complex, typename OrbitalSpinDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename MatsubaraFreqDmn, typename ConcurrencyType>
void computeBareDualGreensFunction(
    const func::function<Complex, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn>>& eps_cg,

    const func::function<Complex, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn>>& eps_tilde,

    const func::function<
        Complex, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>& Delta,

    const func::function<
        Complex, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>& G,
    const ConcurrencyType& concurrency,
    func::function<Complex,
                   func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, MatsubaraFreqDmn>>&

        G0_tilde) {
  // TODO: Add multi-orbital support.
  assert(OrbitalSpinDmn::dmn_size() == 2);

  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, MatsubaraFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency.get_bounds(k_w_dmn_obj);
  int coor[2];

  G0_tilde = 0.;

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, coor);
    const auto k_tilde = coor[0];
    const auto w = coor[1];

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        Complex tmp = -eps_tilde(K1, K2, k_tilde);
        if (K1 == K2)
          tmp += Delta(0, 0, K1, w) + eps_cg(0, 0, K1);

        G0_tilde(K1, K2, k_tilde, w) =
            -G(0, 0, K1, w) * 1. / (G(0, 0, K2, w) + 1. / tmp) * G(0, 0, K2, w);
      }
    }
  }

  concurrency.sum(G0_tilde);
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_COMPUTE_BARE_DUAL_GREENS_FUNCTION_HPP
