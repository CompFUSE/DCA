// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the Gaussian regression method.

#ifndef DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_REGRESSION_HPP
#define DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_REGRESSION_HPP

#include <cmath>
#include <utility>

#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace gp {
// dca::math::gp::

template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
class gaussian_regression {
public:
  typedef typename lhs_dmn_t::parameter_type::scalar_type lhs_scalar_type;
  typedef typename rhs_dmn_t::parameter_type::scalar_type rhs_scalar_type;

  typedef typename lhs_dmn_t::parameter_type::element_type lhs_element_type;
  typedef typename rhs_dmn_t::parameter_type::element_type rhs_element_type;

  gaussian_regression();

  void set_X(linalg::Matrix<scalartype, linalg::CPU>& X_ref) {
    X = X_ref;
  }

  void compute_S(scalartype s_f, scalartype l, scalartype sigma);
  void compute_A(double eps);

  void execute(func::function<scalartype, lhs_dmn_t>& y, func::function<scalartype, rhs_dmn_t>& w);

  void execute(func::function<scalartype, lhs_dmn_t>& y, func::function<scalartype, rhs_dmn_t>& w,
               func::function<scalartype, rhs_dmn_t>& wm, func::function<scalartype, rhs_dmn_t>& wp);

private:
  int Nr;
  int Nc;

  // y = X w
  linalg::Matrix<scalartype, linalg::CPU> X;

  linalg::Matrix<scalartype, linalg::CPU> S;
  linalg::Matrix<scalartype, linalg::CPU> S_inv;

  linalg::Matrix<scalartype, linalg::CPU> A;
  linalg::Matrix<scalartype, linalg::CPU> A_inv;
  linalg::Matrix<scalartype, linalg::CPU> A_inv_X;
};

template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
gaussian_regression<scalartype, lhs_dmn_t, rhs_dmn_t>::gaussian_regression()
    : Nr(rhs_dmn_t::dmn_size()),
      Nc(lhs_dmn_t::dmn_size()),

      X("X", std::pair<int, int>(Nr, Nc)),

      S("S", std::pair<int, int>(Nc, Nc)),
      S_inv("S_inv", std::pair<int, int>(Nc, Nc)),

      A("A", std::pair<int, int>(Nc, Nc)),
      A_inv("A_inv", std::pair<int, int>(Nc, Nc)),
      A_inv_X("A_inv_X", std::pair<int, int>(Nc, Nr)) {}

// Rasmussen and Williams, p. 19, eq. 2.31.
template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
void gaussian_regression<scalartype, lhs_dmn_t, rhs_dmn_t>::compute_S(scalartype s_f, scalartype l,
                                                                      scalartype sigma) {
  for (int i = 0; i < Nc; i++) {
    for (int j = 0; j < Nc; j++) {
      S(i, j) = 0;

      lhs_element_type x_i = lhs_dmn_t::get_elements()[i];
      lhs_element_type x_j = lhs_dmn_t::get_elements()[j];

      S(i, j) += (s_f * s_f) * std::exp(-0.5 * math::util::distance2(x_i, x_j) / (l * l));
    }

    S(i, i) += sigma * sigma;
  }
}

// Rasmussen and Williams, p. 9.
template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
void gaussian_regression<scalartype, lhs_dmn_t, rhs_dmn_t>::compute_A(double eps) {
  linalg::matrixop::pseudoInverse(S, S_inv);

  linalg::Matrix<scalartype, linalg::CPU> Xt_X("Xt_X", std::pair<int, int>(Nc, Nc));

  linalg::matrixop::gemm('C', 'N', X, X, Xt_X);

  for (int j = 0; j < Nc; j++)
    for (int i = 0; i < Nc; i++)
      A(i, j) = Xt_X(i, j) + S(i, j);

  linalg::matrixop::pseudoInverse(A, A_inv);

  linalg::matrixop::gemm('N', 'C', A_inv, X, A_inv_X);
}

template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
void gaussian_regression<scalartype, lhs_dmn_t, rhs_dmn_t>::execute(
    func::function<scalartype, lhs_dmn_t>& y, func::function<scalartype, rhs_dmn_t>& w) {
  w = 0;
  for (int j = 0; j < Nr; j++)
    for (int i = 0; i < Nc; i++)
      w(i) += A_inv_X(i, j) * y(j);
}

template <typename scalartype, typename lhs_dmn_t, typename rhs_dmn_t>
void gaussian_regression<scalartype, lhs_dmn_t, rhs_dmn_t>::execute(
    func::function<scalartype, lhs_dmn_t>& y, func::function<scalartype, rhs_dmn_t>& w,
    func::function<scalartype, rhs_dmn_t>& wm, func::function<scalartype, rhs_dmn_t>& wp) {
  w = 0;
  for (int j = 0; j < Nr; j++)
    for (int i = 0; i < Nc; i++)
      w(i) += A_inv_X(i, j) * y(j);

  for (int i = 0; i < Nc; i++)
    wm(i) = w(i) - A_inv(i, i);

  for (int i = 0; i < Nc; i++)
    wp(i) = w(i) + A_inv(i, i);
}

}  // gp
}  // math
}  // dca

#endif  // DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_REGRESSION_HPP
