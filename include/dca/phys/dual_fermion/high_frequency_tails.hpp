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
#include <iostream>
#include <limits>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, typename OtherDmns>
class HighFrequencyTails {
public:
  using SpFreqType = phys::domains::frequency_domain;
  using SpFreqDmn = func::dmn_0<SpFreqType>;

  using TpFreqType = phys::domains::vertex_frequency_domain<phys::domains::COMPACT>;
  using TpFreqDmn = func::dmn_0<TpFreqType>;

  using TailFreqType = func::ReducedDomain<TpFreqDmn>;
  using TailFreqDmn = func::dmn_0<TailFreqType>;

  HighFrequencyTails(
      const Concurrency& concurrency, const int tail_freqs,
      const func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TpFreqDmn>>& Sigma_tp_freq,
      func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, SpFreqDmn>>& Sigma_sp_freq)
      : concurrency_(concurrency),
        tail_freqs_(tail_freqs),
        Sigma_tp_freq_(Sigma_tp_freq),
        Sigma_sp_freq_(Sigma_sp_freq) {
    // Check tail_freqs parameter.
    if (tail_freqs > TpFreqDmn::dmn_size() / 2)
      throw std::invalid_argument(
          "Number of tail frequencies for fitting cannot be larger than the number of positive "
          "frequencies.");

    if (!TailFreqType::initialized())
      TailFreqType::initialize(0, tail_freqs);
  }

  func::util::Difference compute(bool verbose = false);

private:
  const Concurrency& concurrency_;
  const int tail_freqs_;

  const func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TpFreqDmn>>& Sigma_tp_freq_;
  func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, SpFreqDmn>>& Sigma_sp_freq_;

  func::function<Scalar, OtherDmns> A_;
  func::function<Scalar, OtherDmns> B_;
};

template <typename Scalar, typename Concurrency, typename OtherDmns>
func::util::Difference HighFrequencyTails<Scalar, Concurrency, OtherDmns>::compute(bool verbose) {
  const auto& sp_freqs = SpFreqType::get_elements();
  const auto& tp_freqs = TpFreqType::get_elements();

  // vertex_frequency_domain<COMPACT> cannot be initialized larger than frequency_domain.
  // TODO: Test this branch.
  if (TpFreqDmn::dmn_size() == SpFreqDmn::dmn_size()) {
    // Only copy elements in this case.
    for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
      assert(sp_freqs[w_ind] == tp_freqs[w_ind]);
      for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
        Sigma_sp_freq_(o_ind, w_ind) = Sigma_tp_freq_(o_ind, w_ind);
      }
    }
  }

  // Compute coefficients A and B.
  A_ = 0;
  B_ = 0;

  for (int w_ind = 0; w_ind < tail_freqs_; ++w_ind) {
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      const auto w = tp_freqs[w_ind];
      B_(o_ind) += Sigma_tp_freq_(o_ind, w_ind).real() * w * w;
      B_(o_ind) += Sigma_tp_freq_(o_ind, TpFreqDmn::dmn_size() - 1 - w_ind).real() * w * w;

      A_(o_ind) -= Sigma_tp_freq_(o_ind, w_ind).imag() * w;
      A_(o_ind) += Sigma_tp_freq_(o_ind, TpFreqDmn::dmn_size() - 1 - w_ind).imag() * w;
    }
  }

  A_ /= 2 * tail_freqs_;
  B_ /= 2 * tail_freqs_;

  // Copy tp part.
  for (int w_tp_ind = 0; w_tp_ind < TpFreqDmn::dmn_size(); ++w_tp_ind) {
    const auto w_sp_ind = TpFreqType::get_corresponding_frequency_domain_index()[w_tp_ind];
    assert(std::abs(sp_freqs[w_sp_ind] - tp_freqs[w_tp_ind]) <
           100 * std::numeric_limits<Scalar>::epsilon());

    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      Sigma_sp_freq_(o_ind, w_sp_ind) = Sigma_tp_freq_(o_ind, w_tp_ind);
    }
  }

  // Compute sp part from fit.
  // Negative frequencies.
  for (int w_ind = 0; w_ind < SpFreqDmn::dmn_size() / 2 - TpFreqDmn::dmn_size() / 2; ++w_ind) {
    const auto w = sp_freqs[w_ind];
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      Sigma_sp_freq_(o_ind, w_ind).real(B_(o_ind) / (w * w));
      Sigma_sp_freq_(o_ind, w_ind).imag(-A_(o_ind) / w);
    }
  }

  // Positive frequencies.
  for (int w_ind = SpFreqDmn::dmn_size() / 2 + TpFreqDmn::dmn_size() / 2;
       w_ind < SpFreqDmn::dmn_size(); ++w_ind) {
    const auto w = sp_freqs[w_ind];
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      Sigma_sp_freq_(o_ind, w_ind).real(B_(o_ind) / (w * w));
      Sigma_sp_freq_(o_ind, w_ind).imag(-A_(o_ind) / w);
    }
  }

  // Check quality of fit.
  func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TailFreqDmn>> data;
  func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TailFreqDmn>> fit;

  for (int w_ind = 0; w_ind < tail_freqs_; ++w_ind) {
    for (int o_ind = 0; o_ind < OtherDmns::dmn_size(); ++o_ind) {
      const auto w = tp_freqs[w_ind];
      fit(o_ind, w_ind) = std::complex<Scalar>(B_(o_ind) / (w * w), -A_(o_ind) / w);
      data(o_ind, w_ind) = Sigma_tp_freq_(o_ind, w_ind);
    }
  }
  const auto diff = func::util::difference(data, fit);
  if (verbose) {
    std::cout << "Relative errors of fit of dual-self high-frequency tails:"
              << "\n"
              << "l1    = " << diff.l1 << "\n"
              << "l2    = " << diff.l2 << "\n"
              << "l_inf = " << diff.l_inf << std::endl;
  }

  return diff;
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_HIGH_FREQUENCY_TAILS_HPP
