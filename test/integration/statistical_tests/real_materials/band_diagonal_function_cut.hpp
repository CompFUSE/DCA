// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Description: this file defines the domains of the d.o.f used for the statistic test and the
//              domain of the corresponding covariance matrix.
//              It also provides a method to cut the tested function accordingly.

#ifndef TEST_MATH_STATISTICAL_TESTING_BAND_DIAGONAL_FUNCTION_CUT_HPP
#define TEST_MATH_STATISTICAL_TESTING_BAND_DIAGONAL_FUNCTION_CUT_HPP

#include "dca/function/function.hpp"
#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/reduced_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace testing {
namespace details {
// dca::testing::details::
using Wdmn = func::dmn_0<phys::domains::frequency_domain>;
using Bdmn = func::dmn_0<phys::domains::electron_band_domain>;
using Sdmn = func::dmn_0<phys::domains::electron_spin_domain>;
using Nu = func::dmn_variadic<Bdmn, Sdmn>;  // orbital-spin index

using Rdmn3 = func::dmn_0<phys::domains::cluster_domain<
    double, 3, phys::domains::CLUSTER, phys::domains::REAL_SPACE, phys::domains::BRILLOUIN_ZONE>>;
}  // details
// dca::testing::

using Complex = std::complex<double>;
using DmnWcut = dca::func::ReducedDomain<details::Wdmn>;
using RealImagDmn = func::dmn_0<func::dmn<2, std::string>>;
// RealImagDmn::get_elements() = std::vector<std::string>{"real_part", "imaginary_part"};

template <class SpaceDmn = details::Rdmn3>
using SigmaDomain = func::dmn_variadic<details::Nu, details::Nu, SpaceDmn, details::Wdmn>;
// Stores  n_w frequencies, all k or r points and  bands (no spin).
template <class SpaceDmn = details::Rdmn3>
using SigmaCutDomain = func::dmn_variadic<details::Bdmn, SpaceDmn, func::dmn_0<DmnWcut>, RealImagDmn>;
template <class SpaceDmn = details::Rdmn3>
using CovarianceDomain = func::dmn_variadic<SigmaCutDomain<SpaceDmn>, SigmaCutDomain<SpaceDmn>>;

template <class SpaceDmn = details::Rdmn3>
func::function<double, SigmaCutDomain<SpaceDmn>> cutFrequencyAndOffDiagBand(
    func::function<Complex, SigmaDomain<SpaceDmn>>& f, const int n_w) {
  static int last_n = 0;
  if (last_n != n_w) {
    const int half_point = details::Wdmn::dmn_size() / 2;
    DmnWcut::initialize(half_point, half_point + n_w);
    last_n = n_w;
  }
  const int w0 = details::Wdmn::dmn_size() / 2;
  assert(n_w < w0);
  func::function<double, SigmaCutDomain<SpaceDmn>> result(f.get_name() + "_cut");

  for (int b = 0; b < details::Bdmn::dmn_size(); ++b)
    for (int r = 0; r < SpaceDmn::dmn_size(); ++r)
      for (int iw = 0; iw < n_w; ++iw) {
        result(b, r, iw, 0) = f(b, 0, b, 0, r, w0 + iw).real();
        result(b, r, iw, 1) = f(b, 0, b, 0, r, w0 + iw).imag();
      }

  return result;
}

// r-valued argument version
template <class SpaceDmn>
func::function<double, SigmaCutDomain<SpaceDmn>> cutFrequencyAndOffDiagBand(
    func::function<Complex, SigmaDomain<SpaceDmn>>&& f, const int n_w) {
  return cutFrequencyAndOffDiagBand(f, n_w);
}

}  // testing
}  // dca

#endif  // TEST_MATH_STATISTICAL_TESTING_BAND_DIAGONAL_FUNCTION_CUT_HPP
