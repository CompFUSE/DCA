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

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, typename SpFreqDmn, typename TpFreqDmn, typename OtherDmns>
void highFrequencyTails(
    const Concurrency& concurrency, const int tail_freqs,
    const func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, TpFreqDmn>>& Sigma_tp_freq,
    func::function<std::complex<Scalar>, func::dmn_variadic<OtherDmns, SpFreqDmn>> Sigma_sp_freq,
    Scalar tolerance = 1.e-6) {
  // Check domain sizes.

  // Check tail_freqs parameter.

  // Compute coefficients A and B.

  // Check quality of coefficients.

  // Copy tp part.

  // Compute sp part.
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_HIGH_FREQUENCY_TAILS_HPP
