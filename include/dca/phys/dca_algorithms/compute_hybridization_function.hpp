// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a method to compute the hybridization function,
//     \Delta(\mathbf{K}, \omega_n) = i \omega_n + \mu - \bar{\varepsilon}_{\mathbf{K}}
//                                    - [G_c^0(\mathbf{K}, \omega_n)]^{-1} .

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_HYBDRIDIZATION_FUNCTION_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_HYBDRIDIZATION_FUNCTION_HPP

#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <typename Scalar, typename OrbitalSpinDmn, typename KClusterDmn, typename MatsubaraFreqDmn,
          typename ConcurrencyType>
void computeHybridizationFunction(
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn>>& eps_K_cg,
    const func::function<std::complex<Scalar>, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn,
                                                                  MatsubaraFreqDmn>>& G0_cluster_K_wn,
    const Scalar mu, const ConcurrencyType& concurrency,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KClusterDmn, MatsubaraFreqDmn>>&
        Delta_K_wn) {
  // Work space for the inverse of the bare cluster Green's function.
  linalg::Matrix<std::complex<Scalar>, linalg::CPU> G0_cluster_inv("G0_cluster_inv",
                                                                   OrbitalSpinDmn::dmn_size());
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<std::complex<Scalar>, linalg::CPU> work;

  const std::complex<Scalar> i(0., 1.);

  // Distribute the work amongst the processes.
  const func::dmn_variadic<KClusterDmn, MatsubaraFreqDmn> K_wn_dmn_obj;
  const std::pair<int, int> bounds = concurrency.get_bounds(K_wn_dmn_obj);
  int coor[2];

  Delta_K_wn = 0.;

  for (int l = bounds.first; l < bounds.second; ++l) {
    K_wn_dmn_obj.linind_2_subind(l, coor);
    const auto K_ind = coor[0];
    const auto wn_ind = coor[1];
    const auto wn_val = MatsubaraFreqDmn::get_elements()[wn_ind];

    // Compute G0_cluster^{-1} for fixed k-vector and Matsubara frequency.
    for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        G0_cluster_inv(m, n) = G0_cluster_K_wn(m, n, K_ind, wn_ind);
      }
    }

    linalg::matrixop::inverse(G0_cluster_inv, ipiv, work);

    for (int n = 0; n < OrbitalSpinDmn::dmn_size(); ++n) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        Delta_K_wn(m, n, K_ind, wn_ind) = -eps_K_cg(m, n, K_ind) - G0_cluster_inv(m, n);
        if (m == n)
          Delta_K_wn(m, n, K_ind, wn_ind) += i * wn_val + mu;
      }
    }
  }

  concurrency.sum(Delta_K_wn);
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_HYBDRIDIZATION_FUNCTION_HPP
