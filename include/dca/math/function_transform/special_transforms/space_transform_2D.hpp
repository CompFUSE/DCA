// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class performs the 2D transform from space to momentum used by the tp accumulation.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D

#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <class RDmn, class KDmn, typename Real = double>
class SpaceTransform2D {
private:
  using Complex = std::complex<Real>;

public:
  // Performs the 2D fourier transform from real to momentum space in place.
  // The transform is equivalent to f(k1, k2) = \sum_{r1, r2} exp(i(k1 * r1 - k2 * r2)) f(r1, r2)
  // In/Out: f
  template <class... OtherDmn>
  static void execute(func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmn...>>& f);

  // Specialized version for tp accumulator.
  // Performs the same transformation as above, and the output is rearranged and stored in f_output.
  // In/Out: f_input
  // Out: f_output
  template <class BDmn, class SDmn, class W1Dmn, class W2Dmn>
  static void execute(
      func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, W1Dmn, W2Dmn>>& f_input,
      func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, W1Dmn, W2Dmn>>& f_output);

protected:
  static const linalg::Matrix<Complex, linalg::CPU>& get_T_matrix();
};

template <class RDmn, class KDmn, typename Real>
template <class... OtherDmn>
void SpaceTransform2D<RDmn, KDmn, Real>::execute(
    func::function<Complex, func::dmn_variadic<RDmn, RDmn, OtherDmn...>>& f_input) {
  const int other_size = dca::util::product(OtherDmn::dmn_size()...);
  const int nc = RDmn::dmn_size();
  const int nc2 = nc * nc;
  linalg::Matrix<Complex, linalg::CPU> tmp(nc);
  const auto& T = get_T_matrix();

  for (int l = 0; l < other_size; ++l) {
    linalg::MatrixView<Complex, linalg::CPU> f_r_r(&f_input(l * nc2), nc);
    // f(k1,k2) = \sum exp(i(k1 * r1 - k2 *r2)) f(r1, r2)
    linalg::matrixop::gemm(T, f_r_r, tmp);
    linalg::matrixop::gemm('N', 'C', Complex(1), tmp, T, Complex(0), f_r_r);
  }
}

template <class RDmn, class KDmn, typename Real>
template <class BDmn, class SDmn, class W1Dmn, class W2Dmn>
void SpaceTransform2D<RDmn, KDmn, Real>::execute(
    func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, W1Dmn, W2Dmn>>& f_input,
    func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, W1Dmn, W2Dmn>>& f_output) {
  assert(SDmn::dmn_size() == 2);
  const int nc = RDmn::dmn_size();
  linalg::Matrix<Complex, linalg::CPU> tmp(nc);
  const Complex norm = Complex(1. / nc);
  const auto& T = get_T_matrix();

  for (int w2 = 0; w2 < W2Dmn::dmn_size(); ++w2)
    for (int w1 = 0; w1 < W1Dmn::dmn_size(); ++w1)
      for (int s = 0; s < 2; ++s)
        for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
          for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
            linalg::MatrixView<Complex, linalg::CPU> f_r_r(&f_input(0, 0, b1, b2, s, w1, w2), nc);
            // f(k1,k2) = \sum exp(i(k1 * r1 - k2 *r2)) f(r1, r2)
            linalg::matrixop::gemm(T, f_r_r, tmp);
            linalg::matrixop::gemm('N', 'C', norm, tmp, T, Complex(0), f_r_r);

            for (int k2 = 0; k2 < nc; ++k2)
              for (int k1 = 0; k1 < nc; ++k1)
                f_output(b1, b2, s, k1, k2, w1, w2) = f_r_r(k1, k2);
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

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_TRANSFORM_2D
