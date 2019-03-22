// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class computes the bare dual Green's function.

#ifndef DCA_PHYS_DUAL_FERMION_BARE_DUAL_GREENS_FUNCTION_HPP
#define DCA_PHYS_DUAL_FERMION_BARE_DUAL_GREENS_FUNCTION_HPP

#include <cassert>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Complex, typename Concurrency, typename BandSpinDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename MatsubaraFreqDmn>
class BareDualGreensFunction {
public:
  using Dispersion =
      func::function<Complex, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KClusterDmn>>;
  using NonTransInvariantDispersion =
      func::function<Complex, func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn>>;
  using GF =
      func::function<Complex, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KClusterDmn, MatsubaraFreqDmn>>;
  using DualGF =
      func::function<Complex,
                     func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, MatsubaraFreqDmn>>;

  BareDualGreensFunction(const Concurrency& concurrency, const Dispersion& eps_bar,
                         const NonTransInvariantDispersion& eps_hat, const GF& Delta, const GF& G)
      : concurrency_(concurrency), eps_bar_(eps_bar), eps_hat_(eps_hat), Delta_(Delta), G_(G) {
    // TODO: Multi-orbital support.
    assert(BandSpinDmn::dmn_size() == 2);
  }

  void compute();

  const DualGF& get() const {
    return G0_tilde_;
  }

private:
  const Concurrency& concurrency_;

  const Dispersion& eps_bar_;                   // Coarsegrained dispersion.
  const NonTransInvariantDispersion& eps_hat_;  // Non-translational invariant dispersion.
  const GF& Delta_;                             // Hybridization function.
  const GF& G_;                                 // Cluster Green's function

  DualGF G0_tilde_;  // Bare dual Green's function.
};

template <typename Complex, typename Concurrency, typename BandSpinDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename MatsubaraFreqDmn>
void BareDualGreensFunction<Complex, Concurrency, BandSpinDmn, KClusterDmn, KSuperlatticeDmn,
                            MatsubaraFreqDmn>::compute() {
  // Work space for the inverses.
  linalg::Matrix<Complex, linalg::CPU> m("m", KClusterDmn::dmn_size());
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<Complex, linalg::CPU> work;

  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, MatsubaraFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
  int coor[2];

  G0_tilde_ = 0.;

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, coor);
    const auto k_tilde = coor[0];
    const auto w = coor[1];

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        m(K1, K2) = -eps_hat_(K1, K2, k_tilde);
        if (K1 == K2)
          m(K1, K2) += Delta_(0, 0, K1, w) + eps_bar_(0, 0, K1);
      }
    }

    linalg::matrixop::inverse(m, ipiv, work);

    for (int K = 0; K < KClusterDmn::dmn_size(); ++K) {
      m(K, K) += G_(0, 0, K, w);
    }

    linalg::matrixop::inverse(m, ipiv, work);

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        G0_tilde_(K1, K2, k_tilde, w) = -G_(0, 0, K1, w) * G_(0, 0, K2, w) * m(K1, K2);
      }
    }
  }

  concurrency_.sum(G0_tilde_);
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_BARE_DUAL_GREENS_FUNCTION_HPP
