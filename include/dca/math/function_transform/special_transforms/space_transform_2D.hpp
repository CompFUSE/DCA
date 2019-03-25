// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class implements the general, i.e. for non-translational invariant functions, 2D space Fourier transform.
// See space_to_momentum.hpp for the definition of the corresponding 1D Fourier transform.

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

template <class RDmn, typename Real = double>
class SpaceTransform2D {
protected:
  using Complex = std::complex<Real>;

  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;

  using KDmn = func::dmn_0<typename RDmn::parameter_type::dual_type>;

public:
  // 2D Fourier transform and rearranging of the result.
  // In/Out: f_input
  // Out: f_output
  template <class W1Dmn, class W2Dmn>
  static void execute(
      func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, W1Dmn, W2Dmn>>& f_input,
      func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, W1Dmn, W2Dmn>>& f_output);

  // f(r, r, other-domains) --> f(k, k, other-dmns)
  template <typename OtherDmns>
  static void execute(const func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmns>>& f_in,
                      func::function<Complex, func::dmn_variadic<KDmn, KDmn, OtherDmns>>& f_out);

  // f(k, k, other-domains) --> f(r, r, other-dmns)
  template <typename OtherDmns>
  static void execute(const func::function<Complex, func::dmn_variadic<KDmn, KDmn, OtherDmns>>& f_in,
                      func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmns>>& f_out);

protected:
  static const linalg::Matrix<Complex, linalg::CPU>& get_T_matrix();

  static bool hasPhaseFactors() {
    return SpaceToMomentumTransform<RDmn, KDmn>::hasPhaseFactor();
  }

  static const auto& getPhaseFactors();
};

template <class RDmn, typename Real>
template <class W1Dmn, class W2Dmn>
void SpaceTransform2D<RDmn, Real>::execute(
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

            // f(k1,k2) = \sum exp(i(k1 * r1 - k2 *r2)) f(r1, r2) / Nc
            linalg::matrixop::gemm(T, f_r_r, tmp);
            linalg::matrixop::gemm('N', 'C', norm, tmp, T, Complex(0), f_r_r);

            for (int k2 = 0; k2 < nc; ++k2)
              for (int k1 = 0; k1 < nc; ++k1)
                f_output(b1, b2, s, k1, k2, w1, w2) =
                    f_r_r(k1, k2) * phase_factors(b1, k1) * std::conj(phase_factors(b2, k2));
          }
}

template <class RDmn, typename Real>
template <typename OtherDmns>
void SpaceTransform2D<RDmn, Real>::execute(
    const func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmns>>& f_in,
    func::function<Complex, func::dmn_variadic<KDmn, KDmn, OtherDmns>>& f_out) {
  const int Nc = RDmn::dmn_size();
  const Complex norm = Complex(1. / Nc);

  const auto& T = get_T_matrix();
  linalg::Matrix<Complex, linalg::CPU> tmp(Nc);

  for (int i = 0; i < OtherDmns::dmn_size(); ++i) {
    const auto f_r_r_ptr = linalg::makeConstantView<Complex, linalg::CPU>(&f_in(0, 0, i), Nc);
    linalg::MatrixView<Complex, linalg::CPU> f_k_k(&f_out(0, 0, i), Nc);

    // f(k1,k2) = 1/Nc \sum_{r1, r2} exp(i(k1 * r1 - k2 * r2)) f(r1, r2)
    linalg::matrixop::gemm(T, *f_r_r_ptr, tmp);
    linalg::matrixop::gemm('N', 'C', norm, tmp, T, Complex(0), f_k_k);
  }
}

template <class RDmn, typename Real>
template <typename OtherDmns>
void SpaceTransform2D<RDmn, Real>::execute(
    const func::function<Complex, func::dmn_variadic<KDmn, KDmn, OtherDmns>>& f_in,
    func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmns>>& f_out) {
  const int Nc = RDmn::dmn_size();
  const Complex norm = Complex(1. / Nc);

  const auto& T = get_T_matrix();
  linalg::Matrix<Complex, linalg::CPU> tmp(Nc);

  for (int i = 0; i < OtherDmns::dmn_size(); ++i) {
    const auto f_k_k_ptr = linalg::makeConstantView<Complex, linalg::CPU>(&f_in(0, 0, i), Nc);
    linalg::MatrixView<Complex, linalg::CPU> f_r_r(&f_out(0, 0, i), Nc);

    // f(r1,r2) = 1/Nc \sum_{k1, k2} exp(-i(k1 * r1 - k2 * r2)) f(k1, k2)
    linalg::matrixop::gemm('C', 'N', Complex(1.), T, *f_k_k_ptr, Complex(0.), tmp);
    linalg::matrixop::gemm('N', 'N', norm, tmp, T, Complex(0.), f_r_r);
  }
}

template <class RDmn, typename Real>
const linalg::Matrix<std::complex<Real>, linalg::CPU>& SpaceTransform2D<RDmn, Real>::get_T_matrix() {
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

template <class RDmn, typename Real>
const auto& SpaceTransform2D<RDmn, Real>::getPhaseFactors() {
  static func::function<Complex, func::dmn_variadic<BDmn, KDmn>> phase_factors("Phase factors.");
  static std::once_flag flag;

  // Initialize the phase factors.
  std::call_once(flag, [&]() {
    std::vector<std::vector<Real>> a_vecs;
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
