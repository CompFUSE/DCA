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
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, int dimension>
class DualSelfEnergy {
public:
  using BandDmn = func::dmn_0<phys::domains::electron_band_domain>;

  using FreqExchangeDmn = func::dmn_0<phys::domains::FrequencyExchangeDomain>;
  using TpFreqDmn = func::dmn_0<phys::domains::vertex_frequency_domain<phys::domains::COMPACT>>;
  using DualFreqDmn = func::dmn_0<phys::domains::vertex_frequency_domain<phys::domains::EXTENDED>>;

  using KClusterDmn = typename phys::ClusterDomainAliases<dimension>::KClusterDmn;
  using KSuperlatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;

  using TpGreensFunctionDomain =
      func::dmn_variadic<BandDmn, BandDmn, BandDmn, BandDmn, KClusterDmn, KClusterDmn, KClusterDmn,
                         TpFreqDmn, TpFreqDmn, FreqExchangeDmn>;
  using TpGreensFunction = func::function<std::complex<Scalar>, TpGreensFunctionDomain>;

  using DualGreensFunctionDomain =
      func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, DualFreqDmn>;
  using DualGreensFunction = func::function<std::complex<Scalar>, DualGreensFunctionDomain>;

  DualSelfEnergy(const Concurrency& concurrency, const Scalar beta,
                 const DualGreensFunction& G0_tilde, const TpGreensFunction& Gamma_long_uu,
                 const TpGreensFunction& Gamma_long_ud, const TpGreensFunction& Gamma_tran_ud)
      : concurrency_(concurrency),
        beta_(beta),
        Nc_(KClusterDmn::dmn_size()),
        V_(KSuperlatticeDmn::dmn_size()),
        num_dual_freqs_(DualFreqDmn::dmn_size()),
        num_tp_freqs_(TpFreqDmn::dmn_size()),
        num_exchange_freqs_(FreqExchangeDmn::dmn_size()),
        max_exchange_freq_(FreqExchangeDmn::parameter_type::get_extension_size()),
        G0_tilde_(G0_tilde),
        Gamma_long_uu_(Gamma_long_uu),
        Gamma_long_ud_(Gamma_long_ud),
        Gamma_tran_ud_(Gamma_tran_ud) {
    // TODO: Multi-orbital support.
    assert(BandDmn::dmn_size() == 1);

    // Check consistency of frequency domains.
    assert(num_dual_freqs_ == num_tp_freqs_ + 2 * max_exchange_freq_);
    assert(num_exchange_freqs_ == max_exchange_freq_ + 1);
  }

  // Computes the 1st order contribution.
  void compute1stOrder();

  // Computes the 2nd order contribution.
  void compute2ndOrder() {}

  const DualGreensFunction& get() {
    return Sigma_tilde_;
  }

private:
  const Concurrency& concurrency_;
  const Scalar beta_;

  const int Nc_;
  const int V_;

  const int num_dual_freqs_;
  const int num_tp_freqs_;
  const int num_exchange_freqs_;
  const int max_exchange_freq_;

  // Dual self-energy.
  DualGreensFunction Sigma_tilde_;

  // Bare dual Green's function.
  const DualGreensFunction& G0_tilde_;

  // Two-particle vertex in particle-hole longitudinal up-up channel.
  const TpGreensFunction& Gamma_long_uu_;
  // Two-particle vertex in particle-hole longitudinal up-down channel.
  const TpGreensFunction& Gamma_long_ud_;
  // Two-particle vertex in particle-hole transverse up-down channel.
  const TpGreensFunction& Gamma_tran_ud_;
};

template <typename Scalar, typename Concurrency, int dimension>
void DualSelfEnergy<Scalar, Concurrency, dimension>::compute1stOrder() {
  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, TpFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
  int k_tilde_wn[2];

  Sigma_tilde_ = 0.;

  const Scalar min_1_over_Nc_V_beta = -1. / (Nc_ * V_ * beta_);

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, k_tilde_wn);
    const auto k_tilde = k_tilde_wn[0];
    const auto wn_tp = k_tilde_wn[1];
    const auto wn_sp = wn_tp + max_exchange_freq_;
    assert(DualFreqDmn::get_elements()[wn_sp] == TpFreqDmn::get_elements()[wn_tp]);

    for (int K2 = 0; K2 < Nc_; ++K2) {
      for (int K1 = 0; K1 < Nc_; ++K1) {
        for (int wm_tp = 0; wm_tp < num_tp_freqs_; ++wm_tp) {
          const auto wm_sp = wm_tp + max_exchange_freq_;
          assert(DualFreqDmn::get_elements()[wn_sp] == TpFreqDmn::get_elements()[wn_tp]);

          for (int q_tilde = 0; q_tilde < V_; ++q_tilde) {
            const int k_tilde_plus_q_tilde = KSuperlatticeDmn::parameter_type::add(k_tilde, q_tilde);

            for (int Q = 0; Q < Nc_; ++Q) {
              const int K1_plus_Q = KClusterDmn::parameter_type::add(K1, Q);
              const int K2_plus_Q = KClusterDmn::parameter_type::add(K2, Q);

              Sigma_tilde_(K1, K2, k_tilde, wn_sp) +=
                  min_1_over_Nc_V_beta *
                  (Gamma_long_uu_(0, 0, 0, 0, K1, K2, Q, wn_tp, wm_tp, 0) +
                   Gamma_long_ud_(0, 0, 0, 0, K1, K2, Q, wn_tp, wm_tp, 0)) *
                  G0_tilde_(K1_plus_Q, K2_plus_Q, k_tilde_plus_q_tilde, wm_sp);
            }
          }
        }
      }
    }
  }

  concurrency_.sum(Sigma_tilde_);
}

// template <typename Scalar, typename Concurrency, int dimension>
// void DualSelfEnergy<Scalar, Concurrency, dimension>::compute2ndOrder() {
//   // Distribute the work amongst the processes.
//   const func::dmn_variadic<KSuperlatticeDmn, TpFreqDmn> k_w_dmn_obj;
//   const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
//   int k_tilde_wn[2];

//   Sigma_tilde_ = 0.;

//   const int Nc = KClusterDmn::dmn_size();
//   const int V = KSuperlatticeDmn::dmn_size();

//   const Scalar min_1_over_2_Nc_V_beta_squared = -1. / (2 * Nc * Nc * V * V * beta_ * beta_);

//   for (int l = bounds.first; l < bounds.second; ++l) {
//     k_w_dmn_obj.linind_2_subind(l, k_tilde_wn);
//     const auto k_tilde = k_tilde_wn[0];
//     const auto n = k_tilde_wn[1];

//     for (int K2 = 0; K2 < Nc; ++K2) {
//       for (int K1 = 0; K1 < Nc; ++K1) {
//         // Outer sums.
//         for (int m = 0; m < TpFreqDmn::dmn_size(); ++m) {
//           for (int l = 0; l < FreqExchangeDmn::dmn_size(); ++l) {
//             for (int K2p = 0; K2p < Nc; ++K2p) {
//               for (int K1p = 0; K1p < Nc; ++K1p) {
//                 for (int Q2 = 0; Q2 < Nc; ++Q2) {
//                   const int K2_plus_Q2 = KClusterDmn::parameter_type::add(K2, Q2);
//                   const int K2p_plus_Q2 = KClusterDmn::parameter_type::add(K2p, Q2);

//                   for (int Q1 = 0; Q1 < Nc; ++Q1) {
//                     const int K1_plus_Q1 = KClusterDmn::parameter_type::add(K1, Q1);
//                     const int K1p_plus_Q1 = KClusterDmn::parameter_type::add(K1p, Q1);

//                     const std::complex<Scalar> Gamma_sum_prod =
//                         (Gamma_long_uu_(0, 0, 0, 0, K1, K1p, Q1, n, m, l) *
//                              Gamma_long_uu_(0, 0, 0, 0, K2, K2p, Q2, m, n, l) +
//                          Gamma_long_ud_(0, 0, 0, 0, K1, K1p, Q1, n, m, l) *
//                              Gamma_long_ud_(0, 0, 0, 0, K2, K2p, Q2, m, n, l) +
//                          Gamma_tran_ud_(0, 0, 0, 0, K1, K1p, Q1, n, m, l) *
//                              Gamma_tran_ud_(0, 0, 0, 0, K2, K2p, Q2, m, n, l));

//                     // Inner sums.
//                     for (int k_tilde_p = 0; k_tilde_p < V; ++k_tilde_p) {
//                       for (int q_tilde = 0; q_tilde < V; ++q_tilde) {
//                         const int k_tilde_plus_q_tilde =
//                             KSuperlatticeDmn::parameter_type::add(k_tilde, q_tilde);
//                         const int k_tilde_p_plus_q_tilde =
//                             KSuperlatticeDmn::parameter_type::add(k_tilde_p, q_tilde);

//                         Sigma_tilde_(K1, K2, k_tilde, n) +=
//                             min_1_over_2_Nc_V_beta_squared * Gamma_sum_prod *
//                             G0_tilde_(K1p, K2p, k_tilde_p, m) *
//                             G0_tilde_(K1p_plus_Q1, K2p_plus_Q2, k_tilde_p_plus_q_tilde, m + l) *
//                             G0_tilde_(K1_plus_Q1, K2_plus_Q2, k_tilde_plus_q_tilde, m + l);
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }

//   concurrency_.sum(Sigma_tilde_);
// }

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
