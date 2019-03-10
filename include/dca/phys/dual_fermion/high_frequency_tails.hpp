// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a method to extend the dual self-energy from the smaller two-particle (tp) to
// the larger single-particle (sp) Matsubara freuqency domain.
// The extrapolation is done by fitting the expected high-frequency behavior,
//     \tilde{\Sigma}(i \omega_n >> 0) ~ A/(i \omega_n) + B/\omega_n^2,
// to the tails of the dual self-energy on the tp domain.

#ifndef DCA_PHYS_DUAL_FERMION_HIGH_FREQUENCY_TAILS_HPP
#define DCA_PHYS_DUAL_FERMION_HIGH_FREQUENCY_TAILS_HPP

#include <complex>
#include <limits>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

class HighFrequencyTails {
public:
  using SpFreqType = phys::domains::frequency_domain;
  using SpFreqDmn = func::dmn_0<SpFreqType>;
  using TpFreqType = phys::domains::vertex_frequency_domain<phys::domains::COMPACT>;
  using TpFreqDmn = func::dmn_0<TpFreqType>;

  template <typename Scalar, typename Concurrency, typename OtherDmns>
  static void compute(
      const Concurrency& concurrency, const int tail_freqs,
      const func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TpFreqDmn>>& Sigma_tp_freq,
      func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, SpFreqDmn>>& Sigma_sp_freq,
      Scalar tolerance = 1.e-6) {
    const auto& sp_freqs = SpFreqType::get_elements();
    const auto& tp_freqs = TpFreqType::get_elements();

    // vertex_frequency_domain<COMPACT> cannot be initialized larger than frequency_domain.
    // TODO: Test this branch.
    if (TpFreqDmn::dmn_size() == SpFreqDmn::dmn_size()) {
      // Only copy elements in this case.
      for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
        assert(sp_freqs[w_ind] == tp_freqs[w_ind]);
        for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
          Sigma_sp_freq(o_ind, w_ind) = Sigma_tp_freq(o_ind, w_ind);
        }
      }
    }

    // Check tail_freqs parameter.
    if (tail_freqs > TpFreqDmn::dmn_size() / 2)
      throw std::invalid_argument(
          "Number of tail frequencies for fitting cannot be larger than the number of positive "
          "frequencies.");

    // Compute coefficients A and B.
    func::function<Scalar, OtherDmns> A;
    func::function<Scalar, OtherDmns> B;

    for (int w_ind = 0; w_ind < tail_freqs; ++w_ind) {
      for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
        const auto w = tp_freqs[w_ind];
        B(o_ind) += Sigma_tp_freq(o_ind, w_ind).real() * w * w;
        B(o_ind) += Sigma_tp_freq(o_ind, TpFreqDmn::dmn_size() - 1 - w_ind).real() * w * w;

        A(o_ind) -= Sigma_tp_freq(o_ind, w_ind).imag() * w;
        A(o_ind) += Sigma_tp_freq(o_ind, TpFreqDmn::dmn_size() - 1 - w_ind).imag() * w;
      }
    }

    A /= 2 * tail_freqs;
    B /= 2 * tail_freqs;

    // Copy tp part.
    for (int w_tp_ind = 0; w_tp_ind < TpFreqDmn::dmn_size(); ++w_tp_ind) {
      const auto w_sp_ind = TpFreqType::get_corresponding_frequency_domain_index()[w_tp_ind];
      assert(std::abs(sp_freqs[w_sp_ind] - tp_freqs[w_tp_ind]) <
             100 * std::numeric_limits<Scalar>::epsilon());

      for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
        Sigma_sp_freq(o_ind, w_sp_ind) = Sigma_tp_freq(o_ind, w_tp_ind);
      }
    }

    // Compute sp part from fit.
    // Negative frequencies.
    for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size() / 2 - TpFreqDmn::dmn_size() / 2; ++w_ind) {
      const auto w = sp_freqs[w_ind];
      for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
        Sigma_sp_freq(o_ind, w_ind).real(B(o_ind) / (w * w));
        Sigma_sp_freq(o_ind, w_ind).imag(-A(o_ind) / w);
      }
    }

    // Positive frequencies.
    for (int w_ind = SpFreqDmn::dmn_size() / 2 + TpFreqDmn::dmn_size() / 2;
         w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
      const auto w = sp_freqs[w_ind];
      for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
        Sigma_sp_freq(o_ind, w_ind).real(B(o_ind) / (w * w));
        Sigma_sp_freq(o_ind, w_ind).imag(-A(o_ind) / w);
      }
    }

    // TODO: Check quality of coefficients.
  }
};

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_HIGH_FREQUENCY_TAILS_HPP
