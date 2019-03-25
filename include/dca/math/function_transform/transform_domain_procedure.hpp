// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides transformation characteristics.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_PROCEDURE_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_PROCEDURE_HPP

#include "dca/function/function.hpp"
#include "dca/function/util/real_complex_conversion.hpp"
#include "dca/linalg/linalg.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <int DMN_INDEX>
class TRANSFORM_DOMAIN_PROCEDURE {
public:
  template <typename f_input_t, typename f_output_t>
  static void characterize_transformation(const f_input_t& f_input, const f_output_t& f_output,
                                          int& M, int& K, int& N, int& P);

  template <typename scalartype_1, class domain_input, typename scalartype_2, class domain_output,
            typename scalartype_3>
  static void transform(const func::function<scalartype_1, domain_input>& f_input,
                        func::function<scalartype_2, domain_output>& f_output,
                        const linalg::Matrix<scalartype_3, linalg::CPU>& T);

  template <typename scalartype, class domain_input, class domain_output>
  static void transform(const func::function<scalartype, domain_input>& f_input,
                        func::function<scalartype, domain_output>& f_output,
                        const linalg::Matrix<scalartype, linalg::CPU>& T);

  template <typename scalartype, class domain_input, class domain_output>
  static void transform(const func::function<std::complex<scalartype>, domain_input>& f_input,
                        func::function<std::complex<scalartype>, domain_output>& f_output,
                        const linalg::Matrix<scalartype, linalg::CPU>& T);

  template <typename scalartype, class domain_input, class domain_output>
  static void transform(const func::function<scalartype, domain_input>& f_input,
                        func::function<scalartype, domain_output>& f_output,
                        const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T);

  template <class domain_input, class domain_output>
  static void transform(const func::function<float, domain_input>& f_input,
                        func::function<float, domain_output>& f_output,
                        const linalg::Matrix<double, linalg::CPU>& T);

  template <class domain_input, class domain_output>
  static void transform(const func::function<std::complex<float>, domain_input>& f_input,
                        func::function<std::complex<float>, domain_output>& f_output,
                        const linalg::Matrix<std::complex<double>, linalg::CPU>& T);

  template <typename scalartype, class domain_input, class domain_output>
  static void transform(const func::function<scalartype, domain_output>& f_input,
                        func::function<std::complex<scalartype>, domain_input>& f_output,
                        const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T);

  template <typename scalartype, class domain_input, class domain_output>
  static void transform(const func::function<std::complex<scalartype>, domain_input>& f_input,
                        func::function<scalartype, domain_output>& f_output,
                        const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T);
};

template <int DMN_INDEX>
template <typename f_input_t, typename f_output_t>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::characterize_transformation(const f_input_t& f_input,
                                                                        const f_output_t& f_output,
                                                                        int& M, int& K, int& N,
                                                                        int& P) {
  M = 1;
  for (int l = 0; l < DMN_INDEX; l++)
    M *= f_input[l];

  K = f_input[DMN_INDEX];
  N = f_output[DMN_INDEX];

  P = 1;
  for (int l = DMN_INDEX + 1; l < f_input.signature(); l++)
    P *= f_input[l];
}

template <int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<scalartype, domain_input>& f_input,
    func::function<scalartype, domain_output>& f_output,
    const linalg::Matrix<scalartype, linalg::CPU>& T) {
  int M, K, N, P;
  characterize_transformation(f_input, f_output, M, K, N, P);

  scalartype alpha(1);
  scalartype beta(0);

  if (M == 1) {
    linalg::blas::gemm("N", "N", T.size().first, P, T.size().second, alpha, &T(0, 0),
                       T.leadingDimension(), &f_input(0), f_input[DMN_INDEX], beta, &f_output(0),
                       f_output[DMN_INDEX]);
  }
  else {
    for (int l = 0; l < P; l++) {
      int lin_ind_lhs = M * K * l;
      int lin_ind_rhs = M * N * l;

      linalg::blas::gemm("N", "T", M, N, K, alpha, &f_input(lin_ind_lhs), M, &T(0, 0),
                         T.leadingDimension(), beta, &f_output(lin_ind_rhs), M);
    }
  }
}

template <int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<std::complex<scalartype>, domain_input>& f_input,
    func::function<std::complex<scalartype>, domain_output>& f_output,
    const linalg::Matrix<scalartype, linalg::CPU>& T) {
  int M, K, N, P;
  characterize_transformation(f_input, f_output, M, K, N, P);

  for (int l = 0; l < P; l++) {
    int lin_ind_lhs = M * K * l;
    int lin_ind_rhs = M * N * l;

    scalartype alpha(1);
    scalartype beta(0);

    linalg::blas::gemm("N", "T", 2 * M, N, K, alpha, &real(f_input(lin_ind_lhs)), 2 * M, &T(0, 0),
                       T.leadingDimension(), beta, &real(f_output(lin_ind_rhs)), 2 * M);
  }
}

template <int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<scalartype, domain_input>& f_input,
    func::function<scalartype, domain_output>& f_output,
    const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T) {
  linalg::Matrix<scalartype, linalg::CPU> T_re("T_re", T.size());

  for (int j = 0; j < T.size().second; j++)
    for (int i = 0; i < T.size().first; i++)
      T_re(i, j) = real(T(i, j));

  transform(f_input, f_output, T_re);
}

template <int DMN_INDEX>
template <class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(const func::function<float, domain_input>& f_input,
                                                      func::function<float, domain_output>& f_output,
                                                      const linalg::Matrix<double, linalg::CPU>& T) {
  linalg::Matrix<float, linalg::CPU> T_float("T_re", T.size());

  for (int j = 0; j < T.size().second; j++)
    for (int i = 0; i < T.size().first; i++)
      T_float(i, j) = T(i, j);

  transform(f_input, f_output, T_float);
}

template <int DMN_INDEX>
template <class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<std::complex<float>, domain_input>& f_input,
    func::function<std::complex<float>, domain_output>& f_output,
    const linalg::Matrix<std::complex<double>, linalg::CPU>& T) {
  linalg::Matrix<std::complex<float>, linalg::CPU> T_float("T_re", T.size());

  for (int j = 0; j < T.size().second; j++)
    for (int i = 0; i < T.size().first; i++)
      T_float(i, j) = T(i, j);

  transform(f_input, f_output, T_float);
}

template <int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<scalartype, domain_output>& f_input,
    func::function<std::complex<scalartype>, domain_input>& f_output,
    const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T) {
  const auto f_cmplx = func::util::complex(f_input);

  transform(f_cmplx, f_output, T);
}

template <int DMN_INDEX>
template <typename scalartype, class domain_input, class domain_output>
void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(
    const func::function<std::complex<scalartype>, domain_input>& f_input,
    func::function<scalartype, domain_output>& f_output,
    const linalg::Matrix<std::complex<scalartype>, linalg::CPU>& T) {
  f_output = 0.;

  func::function<scalartype, domain_input> f_in("f_in");
  func::function<scalartype, domain_output> f_out("f_out");

  linalg::Matrix<scalartype, linalg::CPU> T_tmp("T_tmp", T.size());

  {
    for (int i = 0; i < f_input.size(); i++)
      f_in(i) = real(f_input(i));

    for (int j = 0; j < T.size().second; j++)
      for (int i = 0; i < T.size().first; i++)
        T_tmp(i, j) = real(T(i, j));

    transform(f_in, f_out, T_tmp);

    for (int i = 0; i < f_output.size(); i++)
      f_output(i) += f_out(i);
  }

  {
    for (int i = 0; i < f_input.size(); i++)
      f_in(i) = imag(f_input(i));

    for (int j = 0; j < T.size().second; j++)
      for (int i = 0; i < T.size().first; i++)
        T_tmp(i, j) = imag(T(i, j));

    transform(f_in, f_out, T_tmp);

    for (int i = 0; i < f_output.size(); i++)
      f_output(i) -= f_out(i);
  }
}

}  // namespace transform
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_TRANSFORM_DOMAIN_PROCEDURE_HPP
