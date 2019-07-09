// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implements the inverse transformation described in space_to_momentum.hpp.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_MOMENTUM_TO_SPACE_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_MOMENTUM_TO_SPACE_HPP

#include <cassert>
#include <complex>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/domainwise_function_transform.hpp"
#include "dca/math/function_transform/special_transforms/space_to_momentum.hpp"
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
class MomentumToSpaceTransform {
private:
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;

  using Real = typename RDmn::parameter_type::Scalar;
  using Complex = std::complex<Real>;

  using Matrix = linalg::Matrix<Complex, linalg::CPU>;

public:
  // Specialization for complex Sp Green's functions.
  template <class LastDmn>
  static void execute(
      const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_input,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_output);
  static void execute(const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn>>& f_input,
                      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn>>& f_output);

  // Real to complex.
  template <class LastDmn>
  static void execute(
      const func::function<Real, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_input,
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_output) {
    execute(func::util::complex(f_input), f_output);
  }

  // Complex to real.
  // Note: this is not optimized, but it is not performance critical.
  template <class LastDmn>
  static void execute(
      const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_input,
      func::function<Real, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_output) {
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>> f_out_cmplx;
    execute(f_input, f_out_cmplx);
    f_output = std::move(func::util::real(f_out_cmplx, true));
  }

  // Default to old implementation when no band is present.
  // Precondition: DomainInput and DomainOutput don't contain an 'electron_band_domain'.
  template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput,
            class = std::enable_if_t<!dca::util::contained<BDmn, DomainInput>()>>
  static void execute(const func::function<ScalarInp, DomainInput>& f_input,
                      func::function<ScalarOut, DomainOutput>& output);

  static bool hasPhaseFactor();
};

template <class KDmn, class RDmn>
template <class LastDmn>
void MomentumToSpaceTransform<KDmn, RDmn>::execute(
    const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>>& f_input,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>>& f_output) {
  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>> f_input_cpy(f_input);

  if (hasPhaseFactor()) {
    const int n_bands = BDmn::dmn_size();
    const int nc = KDmn::dmn_size();

    auto& phase_factors = SpaceToMomentumTransform<RDmn, KDmn>::getPhaseFactors();

    for (int l = 0; l < LastDmn::dmn_size(); ++l)
      for (int k = 0; k < nc; ++k)
        for (int b2 = 0; b2 < n_bands; ++b2)
          for (int b1 = 0; b1 < n_bands; ++b1)
            for (int s = 0; s < 2; ++s)
              f_input_cpy(b1, s, b2, s, k, l) *= std::conj(phase_factors(b1, b2, k));
  }

  // Execute standard FT
  DomainwiseFunctionTransform<func::dmn_variadic<NuDmn, NuDmn, KDmn, LastDmn>,
                              func::dmn_variadic<NuDmn, NuDmn, RDmn, LastDmn>, typename KDmn::parameter_type,
                              typename RDmn::parameter_type>::execute_on_first(f_input_cpy, f_output);
}

// This implementation is not efficient, but it is anyway only called once on a mall function.
template <class KDmn, class RDmn>
void MomentumToSpaceTransform<KDmn, RDmn>::execute(
    const func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn>>& f_input,
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn>>& f_output) {
  using MockDmn = func::dmn_0<func::dmn<1>>;

  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KDmn, MockDmn>> f_in;
  func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, RDmn, MockDmn>> f_out;
  std::copy_n(f_input.values(), f_input.size(), f_in.values());

  execute(f_in, f_out);

  std::copy_n(f_out.values(), f_input.size(), f_output.values());
}

template <class KDmn, class RDmn>
template <typename ScalarInp, typename ScalarOut, class DomainInput, class DomainOutput, class>
void MomentumToSpaceTransform<KDmn, RDmn>::execute(const func::function<ScalarInp, DomainInput>& f_input,
                                                   func::function<ScalarOut, DomainOutput>& f_output) {
  DomainwiseFunctionTransform<DomainInput, DomainOutput, typename KDmn::parameter_type,
                              typename RDmn::parameter_type>::execute_on_first(f_input, f_output);
}

template <class KDmn, class RDmn>
bool MomentumToSpaceTransform<KDmn, RDmn>::hasPhaseFactor() {
  return SpaceToMomentumTransform<RDmn, KDmn>::hasPhaseFactor();
}

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_MOMENTUM_TO_SPACE_HPP
