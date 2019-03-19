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
#include "dca/math/function_transform/function_transform.hpp"
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
  using ExtFreqDmn = func::dmn_0<phys::domains::vertex_frequency_domain<phys::domains::EXTENDED>>;

  using KClusterDmn = typename phys::ClusterDomainAliases<dimension>::KClusterDmn;
  using KClusterType = typename KClusterDmn::parameter_type;

  using RSuperlatticeDmn = typename phys::ClusterDomainAliases<dimension>::RSpSuperlatticeDmn;
  using KSuperlatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;
  using KSuperlatticeType = typename KSuperlatticeDmn::parameter_type;

  using FTSuperlatticeKtoR =
      dca::math::transform::FunctionTransform<KSuperlatticeDmn, RSuperlatticeDmn>;
  using FTSuperlatticeRtoK =
      dca::math::transform::FunctionTransform<RSuperlatticeDmn, KSuperlatticeDmn>;

  using TpGFDmn = func::dmn_variadic<BandDmn, BandDmn, BandDmn, BandDmn, KClusterDmn, KClusterDmn,
                                     KClusterDmn, TpFreqDmn, TpFreqDmn, FreqExchangeDmn>;
  using TpGF = func::function<std::complex<Scalar>, TpGFDmn>;

  using DualGFExtFreqDmn = func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, ExtFreqDmn>;
  using DualGFExtFreq = func::function<std::complex<Scalar>, DualGFExtFreqDmn>;

  using DualGFTpFreqDmn = func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, TpFreqDmn>;
  using DualGFTpFreq = func::function<std::complex<Scalar>, DualGFTpFreqDmn>;

  DualSelfEnergy(const Concurrency& concurrency, const Scalar beta, const DualGFExtFreq& G0_tilde,
                 const TpGF& Gamma_long_uu, const TpGF& Gamma_long_ud, const TpGF& Gamma_tran_ud)
      : concurrency_(concurrency),

        beta_(beta),
        Nc_(KClusterDmn::dmn_size()),
        V_(KSuperlatticeDmn::dmn_size()),

        K0_(KClusterType::origin_index()),

        max_exchange_freq_(FreqExchangeDmn::parameter_type::get_extension_size()),

        G0_tilde_(G0_tilde),

        Gamma_long_uu_(Gamma_long_uu),
        Gamma_long_ud_(Gamma_long_ud),
        Gamma_tran_ud_(Gamma_tran_ud) {
    // TODO: Multi-orbital support.
    assert(BandDmn::dmn_size() == 1);

    // Check consistency of frequency domains.
    assert(ExtFreqDmn::dmn_size() == TpFreqDmn::dmn_size() + 2 * max_exchange_freq_);
  }

  // Computes the 1st order contribution.
  void compute1stOrder();

  // Computes the 2nd order contribution.
  // Reference implementation.
  void compute2ndOrderReference();
  // Accelerated version using Fourier transformation for the super-lattice variable.
  void compute2ndOrderFT();

  const DualGFTpFreq& get() {
    return Sigma_tilde_;
  }

private:
  inline int minus_w_tp(const int w) const {
    return TpFreqDmn::dmn_size() - 1 - w;
  };

  const Concurrency& concurrency_;
  const Scalar beta_;

  const int Nc_;
  const int V_;

  // Index of origin in momentum space cluster domain.
  const int K0_;

  const int max_exchange_freq_;

  // Dual self-energy.
  DualGFTpFreq Sigma_tilde_;

  // Bare dual Green's function.
  const DualGFExtFreq& G0_tilde_;

  // Two-particle vertex in particle-hole longitudinal up-up channel.
  const TpGF& Gamma_long_uu_;
  // Two-particle vertex in particle-hole longitudinal up-down channel.
  const TpGF& Gamma_long_ud_;
  // Two-particle vertex in particle-hole transverse up-down channel.
  const TpGF& Gamma_tran_ud_;

  // Helper functions for compute2ndOrderFT.
  func::function<std::complex<Scalar>, KSuperlatticeDmn> G1_k;
  func::function<std::complex<Scalar>, KSuperlatticeDmn> G2_k;
  func::function<std::complex<Scalar>, KSuperlatticeDmn> G3_k;

  func::function<std::complex<Scalar>, RSuperlatticeDmn> G1_r;
  func::function<std::complex<Scalar>, RSuperlatticeDmn> G2_r;
  func::function<std::complex<Scalar>, RSuperlatticeDmn> G3_r;
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

    for (int K2 = 0; K2 < Nc_; ++K2) {
      for (int K1 = 0; K1 < Nc_; ++K1) {
        for (int wm_tp = 0; wm_tp < TpFreqDmn::dmn_size(); ++wm_tp) {
          const auto wm_ext = wm_tp + max_exchange_freq_;
          assert(ExtFreqDmn::get_elements()[wm_ext] == TpFreqDmn::get_elements()[wm_tp]);

          for (int q_tilde = 0; q_tilde < V_; ++q_tilde) {
            const int k_tilde_plus_q_tilde = KSuperlatticeType::add(k_tilde, q_tilde);

            for (int Q = 0; Q < Nc_; ++Q) {
              const int K1_plus_Q = KClusterType::add(K1, Q);
              const int K2_plus_Q = KClusterType::add(K2, Q);

              Sigma_tilde_(K1, K2, k_tilde, wn_tp) +=
                  min_1_over_Nc_V_beta *
                  (Gamma_long_uu_(0, 0, 0, 0, K1, K2, Q, wn_tp, wm_tp, 0) +
                   Gamma_long_ud_(0, 0, 0, 0, K1, K2, Q, wn_tp, wm_tp, 0)) *
                  G0_tilde_(K1_plus_Q, K2_plus_Q, k_tilde_plus_q_tilde, wm_ext);
            }
          }
        }
      }
    }
  }

  concurrency_.sum(Sigma_tilde_);
}

template <typename Scalar, typename Concurrency, int dimension>
void DualSelfEnergy<Scalar, Concurrency, dimension>::compute2ndOrderReference() {
  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, TpFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
  int k_tilde_wn[2];

  // Local variable to temporarily store the 2nd order contribution and to do the concurrency sum.
  DualGFTpFreq Sigma_tilde_2nd;

  const Scalar min_1_over_2_Nc_V_beta_squared = -1. / (2. * Nc_ * Nc_ * V_ * V_ * beta_ * beta_);

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, k_tilde_wn);
    const auto k_tilde = k_tilde_wn[0];
    const auto wn_tp = k_tilde_wn[1];

    for (int K2 = 0; K2 < Nc_; ++K2) {
      const int min_K2 = KClusterType::subtract(K2, K0_);

      for (int K1 = 0; K1 < Nc_; ++K1) {
        const int min_K1 = KClusterType::subtract(K1, K0_);

        // Outer sums.
        for (int wm_tp = 0; wm_tp < TpFreqDmn::dmn_size(); ++wm_tp) {
          const auto wm_ext = wm_tp + max_exchange_freq_;
          assert(ExtFreqDmn::get_elements()[wm_ext] == TpFreqDmn::get_elements()[wm_tp]);

          for (int l = -FreqExchangeDmn::dmn_size() + 1; l < FreqExchangeDmn::dmn_size(); ++l) {
            for (int K2p = 0; K2p < Nc_; ++K2p) {
              const int min_K2p = KClusterType::subtract(K2p, K0_);

              for (int K1p = 0; K1p < Nc_; ++K1p) {
                const int min_K1p = KClusterType::subtract(K1p, K0_);

                for (int Q2 = 0; Q2 < Nc_; ++Q2) {
                  const int min_Q2 = KClusterType::subtract(Q2, K0_);
                  const int K2_plus_Q2 = KClusterType::add(K2, Q2);
                  const int K2p_plus_Q2 = KClusterType::add(K2p, Q2);

                  for (int Q1 = 0; Q1 < Nc_; ++Q1) {
                    const int min_Q1 = KClusterType::subtract(Q1, K0_);
                    const int K1_plus_Q1 = KClusterType::add(K1, Q1);
                    const int K1p_plus_Q1 = KClusterType::add(K1p, Q1);

                    std::complex<Scalar> Gamma_sum_prod = 0.;
                    if (l >= 0)
                      Gamma_sum_prod = Gamma_long_uu_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_long_uu_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l) +
                                       Gamma_long_ud_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_long_ud_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l) +
                                       Gamma_tran_ud_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_tran_ud_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l);

                    // Use symmetry of Gamma under conjugation to obtain values at negative bosonic
                    // (exchange) frequencies.
                    else
                      Gamma_sum_prod =
                          std::conj(Gamma_long_uu_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_long_uu_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l)) +
                          std::conj(Gamma_long_ud_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_long_ud_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l)) +
                          std::conj(Gamma_tran_ud_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_tran_ud_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l));

                    // Inner sums (product of \tilde{G}_0's).
                    for (int kp_tilde = 0; kp_tilde < V_; ++kp_tilde) {
                      for (int q_tilde = 0; q_tilde < V_; ++q_tilde) {
                        const int k_tilde_plus_q_tilde = KSuperlatticeType::add(k_tilde, q_tilde);
                        const int kp_tilde_plus_q_tilde = KSuperlatticeType::add(kp_tilde, q_tilde);

                        Sigma_tilde_2nd(K1, K2, k_tilde, wn_tp) +=
                            min_1_over_2_Nc_V_beta_squared * Gamma_sum_prod *
                            G0_tilde_(K1p, K2p, kp_tilde, wm_ext) *
                            G0_tilde_(K1p_plus_Q1, K2p_plus_Q2, kp_tilde_plus_q_tilde, wm_ext + l) *
                            G0_tilde_(K1_plus_Q1, K2_plus_Q2, k_tilde_plus_q_tilde, wm_ext + l);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  concurrency_.sum(Sigma_tilde_2nd);

  Sigma_tilde_ += Sigma_tilde_2nd;
}

template <typename Scalar, typename Concurrency, int dimension>
void DualSelfEnergy<Scalar, Concurrency, dimension>::compute2ndOrderFT() {
  // Distribute the work amongst the processes.
  const func::dmn_variadic<TpFreqDmn> TpFreqDmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(TpFreqDmn_obj);

  // Local variable to temporarily store the 2nd order contribution and to do the concurrency sum.
  DualGFTpFreq Sigma_tilde_2nd;

  // The factor V^2 will be cancelled by two delta-functions from Fourier-transforming the inner
  // sums over \tilde{k'} and \tilde{q}.
  const Scalar min_1_over_2_Nc_beta_squared = -1. / (2. * Nc_ * Nc_ * beta_ * beta_);

  for (int wn_tp = bounds.first; wn_tp < bounds.second; ++wn_tp) {
    for (int K2 = 0; K2 < Nc_; ++K2) {
      const int min_K2 = KClusterType::subtract(K2, K0_);

      for (int K1 = 0; K1 < Nc_; ++K1) {
        const int min_K1 = KClusterType::subtract(K1, K0_);

        // Outer sums.
        for (int wm_tp = 0; wm_tp < TpFreqDmn::dmn_size(); ++wm_tp) {
          const auto wm_ext = wm_tp + max_exchange_freq_;
          assert(ExtFreqDmn::get_elements()[wm_ext] == TpFreqDmn::get_elements()[wm_tp]);

          for (int l = -FreqExchangeDmn::dmn_size() + 1; l < FreqExchangeDmn::dmn_size(); ++l) {
            for (int K2p = 0; K2p < Nc_; ++K2p) {
              const int min_K2p = KClusterType::subtract(K2p, K0_);

              for (int K1p = 0; K1p < Nc_; ++K1p) {
                const int min_K1p = KClusterType::subtract(K1p, K0_);

                for (int Q2 = 0; Q2 < Nc_; ++Q2) {
                  const int min_Q2 = KClusterType::subtract(Q2, K0_);
                  const int K2_plus_Q2 = KClusterType::add(K2, Q2);
                  const int K2p_plus_Q2 = KClusterType::add(K2p, Q2);

                  for (int Q1 = 0; Q1 < Nc_; ++Q1) {
                    const int min_Q1 = KClusterType::subtract(Q1, K0_);
                    const int K1_plus_Q1 = KClusterType::add(K1, Q1);
                    const int K1p_plus_Q1 = KClusterType::add(K1p, Q1);

                    std::complex<Scalar> Gamma_sum_prod = 0.;
                    if (l >= 0)
                      Gamma_sum_prod = Gamma_long_uu_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_long_uu_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l) +
                                       Gamma_long_ud_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_long_ud_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l) +
                                       Gamma_tran_ud_(0, 0, 0, 0, K1, K1p, Q1, wn_tp, wm_tp, l) *
                                           Gamma_tran_ud_(0, 0, 0, 0, K2, K2p, Q2, wm_tp, wn_tp, l);

                    // Use symmetry of Gamma under conjugation to obtain values at negative bosonic
                    // (exchange) frequencies.
                    else
                      Gamma_sum_prod =
                          std::conj(Gamma_long_uu_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_long_uu_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l)) +
                          std::conj(Gamma_long_ud_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_long_ud_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l)) +
                          std::conj(Gamma_tran_ud_(0, 0, 0, 0, min_K1, min_K1p, min_Q1,
                                                   minus_w_tp(wn_tp), minus_w_tp(wm_tp), -l)) *
                              std::conj(Gamma_tran_ud_(0, 0, 0, 0, min_K2, min_K2p, min_Q2,
                                                       minus_w_tp(wm_tp), minus_w_tp(wn_tp), -l));

                    // The inner sums over \tilde{k'} and \tilde{q} are replaced by super-lattice
                    // Fourier transforms and the loop to update \tilde{k} is moved inside the outer
                    // sums.

                    // Set G1_k, G2_k, G3_k.
                    for (int k_tilde = 0; k_tilde < V_; ++k_tilde) {
                      G1_k(k_tilde) = G0_tilde_(K1p, K2p, k_tilde, wm_ext);
                      G2_k(k_tilde) = G0_tilde_(K1p_plus_Q1, K2p_plus_Q2, k_tilde, wm_ext + l);
                      G3_k(k_tilde) = G0_tilde_(K1_plus_Q1, K2_plus_Q2, k_tilde, wm_ext + l);
                    }

                    // Super-lattice FTs of G1, G2 and G3.
                    FTSuperlatticeKtoR::execute(G1_k, G1_r);
                    FTSuperlatticeKtoR::execute(G2_k, G2_r);
                    FTSuperlatticeKtoR::execute(G3_k, G3_r);

                    // Element-wise multiplication of G1, G2 and G3. Result is stored in G1.
                    for (int r_tilde = 0; r_tilde < V_; ++r_tilde) {
                      G1_r(r_tilde) *= G2_r(r_tilde) * G3_r(r_tilde);
                    }

                    // Inverse super-lattice FT of G1.
                    FTSuperlatticeRtoK::execute(G1_r, G1_k);

                    // Update Sigma_tilde_2nd.
                    for (int k_tilde = 0; k_tilde < V_; ++k_tilde) {
                      Sigma_tilde_2nd(K1, K2, k_tilde, wn_tp) +=
                          min_1_over_2_Nc_beta_squared * Gamma_sum_prod * G1_k(k_tilde);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  concurrency_.sum(Sigma_tilde_2nd);

  Sigma_tilde_ += Sigma_tilde_2nd;
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
