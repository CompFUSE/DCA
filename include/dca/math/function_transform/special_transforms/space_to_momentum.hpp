// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Specialization for momentum to space transforms with phase factors.
// This class implements the transform of single particle functions equivalent to the transformation
// of  fermionic operators
// c(b, k) = \sum_e Exp[i k (r + a[b])] c(b, k) / sqrt[Nc],
// where a[b] describes the position of the orbital b.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TO_MOMENTUM_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TO_MOMENTUM_HPP

#include <cassert>
#include <complex>
#include <mutex>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/domainwise_function_transform.hpp"
#include "dca/linalg/blas/blas3.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <class KDmn, class RDmn>
class MomentumToSpaceTransform;

template <class RDmn, class KDmn>
class SpaceToMomentumTransform {
private:
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;

  using Real = typename RDmn::parameter_type::Scalar;
  using Complex = std::complex<Real>;

public:
  // Specialization for complex Sp Green's functions.
  template <class LastDmn>
  static void execute(
      const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_input,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_output);

  // Real to complex.
  template <class LastDmn>
  static void execute(
      const func::function<Real, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_input,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_output) {
    execute(func::util::complex(f_input), f_output);
  }

  // Complex to real.
  // Note: this is not optimized, but it is not performance critical.
  template <class LastDmn>
  static void execute(
      const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_input,
      func::function<Real, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_output) {
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>> f_out_cmplx;
    execute(f_input, f_out_cmplx);
    f_output = std::move(func::util::real(f_out_cmplx, true));
  }

  // A non specialized transform with an 'electron_band_domain' is forbidden.
  template <typename ScalarA, typename ScalarB, class... PackIn, class OutDmn>
  static void execute(const func::function<ScalarA, func::dmn_variadic<BDmn, PackIn...>>& /*f_in*/,
                      func::function<ScalarB, OutDmn>& /*f_out*/) = delete;
  template <typename ScalarA, typename ScalarB, class... PackIn, class OutDmn>
  static void execute(const func::function<ScalarA, func::dmn_variadic<NuDmn, PackIn...>>& /*f_in*/,
                      func::function<ScalarB, OutDmn>& /*f_out*/) = delete;

  // Default to old implementation when no 'electron_band_domain' is present.
  template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput>
  static void execute(const func::function<ScalarInp, DomainInput>& f_input,
                      func::function<ScalarOut, DomainOutput>& output);

  static bool hasPhaseFactor();

private:
  using Matrix = linalg::Matrix<Complex, linalg::CPU>;
  using PhaseFactors = func::function<Complex, func::dmn_variadic<BDmn, BDmn, KDmn>>;

  static const PhaseFactors& getPhaseFactors();

  friend class MomentumToSpaceTransform<KDmn, RDmn>;
};

template <class RDmn, class KDmn>
template <class LastDmn>
void SpaceToMomentumTransform<RDmn, KDmn>::execute(
    const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_input,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_output) {
  // Execute standard FT
  DomainwiseFunctionTransform<func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>,
                              func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>, typename RDmn::parameter_type,
                              typename KDmn::parameter_type>::execute_on_first(f_input, f_output);

  if (hasPhaseFactor()) {
    const int n_bands = BDmn::dmn_size();
    const int nc = KDmn::dmn_size();

    const auto& phase_factors = getPhaseFactors();

    for (int l = 0; l < LastDmn::dmn_size(); ++l) {
      for (int k = 0; k < nc; ++k)
        for (int b2 = 0; b2 < n_bands; ++b2)
          for (int b1 = 0; b1 < n_bands; ++b1)
            for (int s = 0; s < 2; ++s)
              f_output(b1, s, b2, s, k, l) *= phase_factors(b1, b2, k);
    }
  }
}

template <class RDmn, class KDmn>
template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput>
void SpaceToMomentumTransform<RDmn, KDmn>::execute(const func::function<ScalarInp, DomainInput>& f_input,
                                                   func::function<ScalarOut, DomainOutput>& f_output) {
  DomainwiseFunctionTransform<DomainInput, DomainOutput, typename RDmn::parameter_type,
                              typename KDmn::parameter_type>::execute_on_first(f_input, f_output);
}

template <class RDmn, class KDmn>
bool SpaceToMomentumTransform<RDmn, KDmn>::hasPhaseFactor() {
  auto check_a_vecs = [&]() {
    for (const auto& dmn_elem : BDmn::parameter_type::get_elements())
      for (const auto a : dmn_elem.a_vec)
        if (a != 0)
          return true;
    return false;
  };

  const static bool result = check_a_vecs();
  return result;
}

template <class RDmn, class KDmn>
const typename SpaceToMomentumTransform<RDmn, KDmn>::PhaseFactors& SpaceToMomentumTransform<
    RDmn, KDmn>::getPhaseFactors() {
  static PhaseFactors phase_factors("Phase factors.");
  static std::once_flag flag;

  // Initialize the phase factors.
  std::call_once(flag, [&]() {
    std::vector<std::vector<Real>> a_vecs;
    for (const auto& elem : BDmn::get_elements())
      a_vecs.push_back(elem.a_vec);

    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      const auto& k_vec = KDmn::get_elements()[k];
      for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
        for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1)
          phase_factors(b1, b2, k) = std::exp(
              Complex(0., util::innerProduct(k_vec, util::subtract(a_vecs[b2], a_vecs[b1]))));
    }
  });

  return phase_factors;
}

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TO_MOMENTUM_HPP
