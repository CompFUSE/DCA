// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class performs the 2D transform from space to momentum used by the tp accumulation on
// non translational invariant Greens function. See space_to_momentum.hpp for a definition of the
// transformation.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D

#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/math/function_transform/special_transforms/space_to_momentum.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <class RDmn, class KDmn, typename Real = double>
class SpaceTransform2D {
protected:
  using Complex = std::complex<Real>;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;

public:
  // Apply the 2D Fourier transform defined as
  // f(k1, b1, k2, b2) = \sum_{r1, r2} Exp[i (k1 (r1 + a[b1]) - k2 (r2 + a[b2])] f(r1, b1, r2, b2) / Nc,
  // where a[b] is the displacement vector associated with each band,
  // and rearrange the output domains.
  // In/Out: f_input. The input is overwritten with a partial result.
  // Out: f_output
  template <class W1Dmn, class W2Dmn>
  static void execute(
      func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, W1Dmn, W2Dmn>>& f_input,
      func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, W1Dmn, W2Dmn>>& f_output);

protected:
  static const linalg::Matrix<Complex, linalg::CPU>& get_T_matrix();

  static bool hasPhaseFactors() {
    return SpaceToMomentumTransform<RDmn, KDmn>::hasPhaseFactor();
  }

  static const auto& getPhaseFactors();
};

template <class RDmn, class KDmn, typename Real>
template <class W1Dmn, class W2Dmn>
void SpaceTransform2D<RDmn, KDmn, Real>::execute(
    func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, W1Dmn, W2Dmn>>& f_input,
    func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, W1Dmn, W2Dmn>>& f_output) {
  assert(SDmn::dmn_size() == 2);
  const int nc = RDmn::dmn_size();
  linalg::Matrix<Complex, linalg::CPU> tmp(nc);
  const Complex norm = Complex(1. / nc);
  const auto& T = get_T_matrix();
  const auto& phase_factors = getPhaseFactors();

  for (int w2 = 0; w2 < W2Dmn::dmn_size(); ++w2)
    for (int w1 = 0; w1 < W1Dmn::dmn_size(); ++w1)
      for (int s = 0; s < 2; ++s)
        for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
          for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
            linalg::MatrixView<Complex, linalg::CPU> f_r_r(&f_input(0, 0, b1, b2, s, w1, w2), nc);

            // f(k1,k2) = \sum_{r1, r2} exp(i(k1 * r1 - k2 *r2)) f(r1, r2) / Nc
            linalg::matrixop::gemm(T, f_r_r, tmp);
            linalg::matrixop::gemm('N', 'C', norm, tmp, T, Complex(0), f_r_r);

            // f(k1, k2) *= Exp[i (k1 a[b1] - k2 a[b2])]
            for (int k2 = 0; k2 < nc; ++k2)
              for (int k1 = 0; k1 < nc; ++k1)
                f_output(b1, b2, s, k1, k2, w1, w2) =
                    f_r_r(k1, k2) * phase_factors(b1, k1) * std::conj(phase_factors(b2, k2));
          }
}

template <class RDmn, class KDmn, typename Real>
const linalg::Matrix<std::complex<Real>, linalg::CPU>& SpaceTransform2D<RDmn, KDmn, Real>::get_T_matrix() {
  auto initialize_T_matrix = []() {
    assert(RDmn::dmn_size() == KDmn::dmn_size());
    linalg::Matrix<Complex, linalg::CPU> T(RDmn::dmn_size());
    for (int j = 0; j < RDmn::dmn_size(); ++j) {
      const auto& r = RDmn::parameter_type::get_elements()[j];
      for (int i = 0; i < KDmn::dmn_size(); ++i) {
        const auto& k = KDmn::parameter_type::get_elements()[i];
        T(i, j) = std::exp(Complex(0, util::innerProduct(k, r)));
      }
    }
    return T;
  };

  static const auto T = initialize_T_matrix();
  return T;
}

template <class RDmn, class KDmn, typename Real>
const auto& SpaceTransform2D<RDmn, KDmn, Real>::getPhaseFactors() {
  static func::function<Complex, func::dmn_variadic<BDmn, KDmn>> phase_factors("Phase factors.");
  static std::once_flag flag;

  // Initialize the phase factors Exp[i k a[b]].
  std::call_once(flag, [&]() {
    std::vector<std::vector<double>> a_vecs;
    for (const auto& elem : BDmn::get_elements())
      a_vecs.push_back(elem.a_vec);

    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      const auto& k_vec = KDmn::get_elements()[k];
      for (int b = 0; b < BDmn::dmn_size(); ++b)
        phase_factors(b, k) = std::exp(Complex(0., util::innerProduct(k_vec, a_vecs[b])));
    }
  });

  return phase_factors;
}

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D
