// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class implements the dual self-energy up to second order.

#ifndef DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
#define DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP

#include <cassert>
#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, typename BandDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename TpFreqDmn, typename FreqExchangeDmn>
class DualSelfEnergy {
public:
  using TpGreensFunctionDomain =
      func::dmn_variadic<BandDmn, BandDmn, BandDmn, BandDmn, KClusterDmn, KClusterDmn, KClusterDmn,
                         TpFreqDmn, TpFreqDmn, FreqExchangeDmn>;
  using TpGreensFunction = func::function<std::complex<Scalar>, TpGreensFunctionDomain>;

  using DualGreensFunctionDomain =
      func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, TpFreqDmn>;
  using DualGreensFunction = func::function<std::complex<Scalar>, DualGreensFunctionDomain>;

  DualSelfEnergy(const Concurrency& concurrency, const Scalar beta,
                 const DualGreensFunction& G0_tilde, const TpGreensFunction& Gamma_long_uu,
                 const TpGreensFunction& Gamma_long_ud, const TpGreensFunction& Gamma_trans_ud)
      : concurrency_(concurrency),
        beta_(beta),
        G0_tilde_(G0_tilde),
        Gamma_long_uu_(Gamma_long_uu),
        Gamma_long_ud_(Gamma_long_ud),
        Gamma_trans_ud_(Gamma_trans_ud) {
    // TODO: Multi-orbital support.
    assert(BandDmn::dmn_size() == 1);
  }

  // Computes the 1st order contribution.
  void compute1stOrder();

  // Computes the 2nd order contribution.
  void compute2ndOrder(){};

  const DualGreensFunction& get() {
    return Sigma_tilde_;
  }

private:
  const Concurrency& concurrency_;
  const Scalar beta_;

  // Dual self-energy.
  DualGreensFunction Sigma_tilde_;

  // Bare dual Green's function.
  const DualGreensFunction& G0_tilde_;

  // Two-particle vertex in particle-hole longitudinal up-up channel.
  const TpGreensFunction& Gamma_long_uu_;
  // Two-particle vertex in particle-hole longitudinal up-down channel.
  const TpGreensFunction& Gamma_long_ud_;
  // Two-particle vertex in particle-hole transverse up-down channel.
  const TpGreensFunction& Gamma_trans_ud_;
};

template <typename Scalar, typename Concurrency, typename BandDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename TpFreqDmn, typename FreqExchangeDmn>
void DualSelfEnergy<Scalar, Concurrency, BandDmn, KClusterDmn, KSuperlatticeDmn, TpFreqDmn,
                    FreqExchangeDmn>::compute1stOrder() {
  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, TpFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
  int k_tilde_wn[2];

  Sigma_tilde_ = 0.;

  const double min_1_over_Nc_V_beta =
      -1. / (KClusterDmn::dmn_size() * KSuperlatticeDmn::dmn_size() * beta_);

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, k_tilde_wn);
    const auto k_tilde = k_tilde_wn[0];
    const auto wn = k_tilde_wn[1];

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        for (int wm = 0; wm < TpFreqDmn::dmn_size(); ++wm) {
          for (int q_tilde = 0; q_tilde < KSuperlatticeDmn::dmn_size(); ++q_tilde) {
            const int k_tilde_plus_q_tilde = KSuperlatticeDmn::parameter_type::add(k_tilde, q_tilde);

            for (int Q = 0; Q < KClusterDmn::dmn_size(); ++Q) {
              const int K1_plus_Q = KClusterDmn::parameter_type::add(K1, Q);
              const int K2_plus_Q = KClusterDmn::parameter_type::add(K2, Q);

              Sigma_tilde_(K1, K2, k_tilde, wn) +=
                  min_1_over_Nc_V_beta *
                  (Gamma_long_uu_(0, 0, 0, 0, K1, K2, Q, wn, wm, 0) +
                   Gamma_long_ud_(0, 0, 0, 0, K1, K2, Q, wn, wm, 0)) *
                  G0_tilde_(K1_plus_Q, K2_plus_Q, k_tilde_plus_q_tilde, wm);
            }
          }
        }
      }
    }
  }

  concurrency_.sum(Sigma_tilde_);
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
