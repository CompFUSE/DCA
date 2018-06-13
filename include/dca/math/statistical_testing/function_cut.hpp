// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Description: this file defines the domains of the d.o.f used by the statistical test and the
//              domain of the corresponding covariance matrix.
//              It also provides  methods to cut the tested function accordingly.

#ifndef DCA_MATH_STATISTICAL_TESTING_FUNCTION_CUT_HPP
#define DCA_MATH_STATISTICAL_TESTING_FUNCTION_CUT_HPP

#include <cassert>
#include <string>
#include <vector>

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/reduced_domain.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace math {
namespace util {
namespace details {
// dca::math::util::details::
using Wdmn = func::dmn_0<phys::domains::frequency_domain>;
using Bdmn = func::dmn_0<phys::domains::electron_band_domain>;
using Sdmn = func::dmn_0<phys::domains::electron_spin_domain>;
using Nu = func::dmn_variadic<Bdmn, Sdmn>;
// Default momentum space.
using Kdmn = func::dmn_0<phys::domains::cluster_domain<
    double, 2, phys::domains::CLUSTER, phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;
}  // details
// dca::math::util::

using Complex = std::complex<double>;
using DmnWcut = dca::func::ReducedDomain<details::Wdmn>;
using RealImagDmn = func::dmn_0<func::dmn<2, std::string>>;
template <class SpaceDmn = details::Kdmn>
using SigmaDomain = func::dmn_variadic<details::Nu, details::Nu, SpaceDmn, details::Wdmn>;

// Stores  n_w frequencies, all k or r points and bands (no spin).
template <class SpaceDmn = details::Kdmn>
using SigmaCutDomain =
    func::dmn_variadic<details::Bdmn, details::Bdmn, SpaceDmn, func::dmn_0<DmnWcut>, RealImagDmn>;
template <class SpaceDmn = details::Kdmn>
using CovarianceDomain = func::dmn_variadic<SigmaCutDomain<SpaceDmn>, SigmaCutDomain<SpaceDmn>>;

// Prepares a function in the domain of the Self Energy to be statistically tested.
// It removes the spin sector, unroll the real and imaginary part, and cut keeps only n_w positive
// frequencies.
// In: f, n_w.
// Returns: the cut function.
template <class SpaceDmn = details::Kdmn>
func::function<double, SigmaCutDomain<SpaceDmn>> cutFrequency(
    const func::function<Complex, SigmaDomain<SpaceDmn>>& f, const int n_w) {
  assert(n_w > 0);
  static int last_n = -1;
  if (last_n != n_w) {
    const int half_point = details::Wdmn::dmn_size() / 2;
    DmnWcut::initialize(half_point, half_point + n_w);
    RealImagDmn::parameter_type::set_elements(
        std::vector<std::string>{"real_part", "imaginary_part"});
    last_n = n_w;
  }
  const int w0 = details::Wdmn::dmn_size() / 2;
  assert(n_w < w0);
  func::function<double, SigmaCutDomain<SpaceDmn>> result(f.get_name() + "_cut");

  for (int b1 = 0; b1 < details::Bdmn::dmn_size(); ++b1)
    for (int b2 = 0; b2 < details::Bdmn::dmn_size(); ++b2)
      for (int ik = 0; ik < SpaceDmn::dmn_size(); ++ik)
        for (int iw = 0; iw < n_w; ++iw) {
          result(b1, b2, ik, iw, 0) = f(b1, 0, b2, 0, ik, w0 + iw).real();
          result(b1, b2, ik, iw, 1) = f(b1, 0, b2, 0, ik, w0 + iw).imag();
        }

  return result;
}

// Stores only the diagonal band component of a function whose first two
// indices belong to the band domain.
// In: f.
// Returns: the cut function.
template <typename Scalar, class... Domains>
func::function<Scalar, func::dmn_variadic<details::Bdmn, Domains...>> bandDiagonal(
    const func::function<Scalar, func::dmn_variadic<details::Bdmn, details::Bdmn, Domains...>>& f) {
  func::function<Scalar, func::dmn_variadic<details::Bdmn, Domains...>> result;
  const int n_b = details::Bdmn::dmn_size();
  for (int base_idx_f = 0, base_idx_res = 0; base_idx_f < f.size();
       base_idx_f += n_b * n_b, base_idx_res += n_b)
    for (int b = 0; b < n_b; ++b)
      result(b + base_idx_res) = f(b * (n_b + 1) + base_idx_f);
  return result;
}

// r-value reference argument versions.
template <class SpaceDmn>
func::function<double, SigmaCutDomain<SpaceDmn>> cutFrequency(
    func::function<Complex, SigmaDomain<SpaceDmn>>&& f, const int n_w) {
  return cutFrequency(f, n_w);
}
template <typename Scalar, class... Domains>
func::function<Scalar, func::dmn_variadic<details::Bdmn, Domains...>> bandDiagonal(
    func::function<Scalar, func::dmn_variadic<details::Bdmn, details::Bdmn, Domains...>>&& f) {
  return bandDiagonal(f);
}

}  // util
}  // math
}  // dca

#endif  // DCA_MATH_STATISTICAL_TESTING_FUNCTION_CUT_HPP
