// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the analytical transformation used in the interpolation of the self-energy.

/*
 *  \f{eqnarray*}{
 *     T     [\Sigma] &=& [ \Sigma - \alpha I ]^{-1}
 *     T^{-1}[\Sigma] &=&  T[\Sigma]^{-1} + \alpha I
 *  \f}
 */

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_TRANSFORM_TO_ALPHA_HPP
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_TRANSFORM_TO_ALPHA_HPP

#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "comp_library/linalg/linalg.hpp"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {
class transform_to_alpha {
public:
  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

public:
  template <typename scalar_type, typename k_dmn_t>
  void forward(scalar_type alpha, func::function<std::complex<scalar_type>, k_dmn_t>& f_k,
               func::function<std::complex<scalar_type>, k_dmn_t>& alpha_k);

  template <typename scalar_type, typename k_dmn_t>
  void backward(scalar_type alpha, func::function<std::complex<scalar_type>, k_dmn_t>& f_k,
                func::function<std::complex<scalar_type>, k_dmn_t>& alpha_k);

  template <typename scalar_type, typename k_dmn_t>
  static void forward(
      scalar_type alpha,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& alpha_k_w);

  template <typename scalar_type, typename k_dmn_t>
  static void backward(
      scalar_type alpha,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& alpha_k_w);

  template <typename scalar_type, typename k_dmn_t, typename w_dmn_t>
  static void forward(
      scalar_type alpha,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& alpha_k_w);

  template <typename scalar_type, typename k_dmn_t, typename w_dmn_t>
  static void backward(
      scalar_type alpha,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_k_w,
      func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& alpha_k_w);
};

template <typename scalar_type, typename k_dmn_t>
void transform_to_alpha::forward(scalar_type alpha,
                                 func::function<std::complex<scalar_type>, k_dmn_t>& f_k,
                                 func::function<std::complex<scalar_type>, k_dmn_t>& alpha_k) {
  std::complex<scalar_type> I(0., alpha);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
    alpha_k(k_ind) = 1. / (f_k(k_ind) - I);
}

template <typename scalar_type, typename k_dmn_t>
void transform_to_alpha::backward(scalar_type alpha,
                                  func::function<std::complex<scalar_type>, k_dmn_t>& f_k,
                                  func::function<std::complex<scalar_type>, k_dmn_t>& alpha_k) {
  std::complex<scalar_type> I(0., alpha);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
    f_k(k_ind) = 1. / alpha_k(k_ind) + I;
}

template <typename scalar_type, typename k_dmn_t>
void transform_to_alpha::forward(
    scalar_type alpha,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& alpha_k_w) {
  std::complex<scalar_type> I(0., alpha);

  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> f_matrix(
      "f_matrix", std::pair<int, int>(N, N));
  LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(f_matrix);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        f_matrix(i, j) = f_k_w(i, j, k_ind);

    for (int i = 0; i < N; i++)
      f_matrix(i, i) -= I;

    // LIN_ALG::GEINV<dca::linalg::CPU>::execute_on_Green_function_matrix(f_matrix);
    geinv_obj.execute(f_matrix);

    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        alpha_k_w(i, j, k_ind) = f_matrix(i, j);
  }
}

template <typename scalar_type, typename k_dmn_t>
void transform_to_alpha::backward(
    scalar_type alpha,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& f_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t>>& alpha_k_w) {
  std::complex<scalar_type> I(0., alpha);

  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> f_matrix(
      "f_matrix", std::pair<int, int>(N, N));
  LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(f_matrix);

  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        f_matrix(i, j) = alpha_k_w(i, j, k_ind);

    // LIN_ALG::GEINV<dca::linalg::CPU>::execute_on_Green_function_matrix(f_matrix);
    geinv_obj.execute(f_matrix);

    for (int i = 0; i < N; i++)
      f_matrix(i, i) += I;

    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        f_k_w(i, j, k_ind) = f_matrix(i, j);
  }
}

template <typename scalar_type, typename k_dmn_t, typename w_dmn_t>
void transform_to_alpha::forward(
    scalar_type alpha,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& alpha_k_w) {
  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> f_matrix(
      "f_matrix", std::pair<int, int>(N, N));
  LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(f_matrix);

  for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
    scalar_type factor = w_dmn_t::get_elements()[w_ind] > 0 ? 1 : -1;

    std::complex<scalar_type> I(0., factor * alpha);

    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
          f_matrix(i, j) = f_k_w(i, j, k_ind, w_ind);

      for (int i = 0; i < N; i++)
        f_matrix(i, i) -= I;

      // LIN_ALG::GEINV<dca::linalg::CPU>::execute_on_Green_function_matrix(f_matrix);
      geinv_obj.execute(f_matrix);

      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
          alpha_k_w(i, j, k_ind, w_ind) = f_matrix(i, j);
    }
  }
}

template <typename scalar_type, typename k_dmn_t, typename w_dmn_t>
void transform_to_alpha::backward(
    scalar_type alpha,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& f_k_w,
    func::function<std::complex<scalar_type>, func::dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>& alpha_k_w) {
  int N = nu::dmn_size();

  dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU> f_matrix(
      "f_matrix", std::pair<int, int>(N, N));
  LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<scalar_type>> geinv_obj(f_matrix);

  for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++) {
    scalar_type factor = w_dmn_t::get_elements()[w_ind] > 0 ? 1 : -1;

    std::complex<scalar_type> I(0., factor * alpha);

    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++) {
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
          f_matrix(i, j) = alpha_k_w(i, j, k_ind, w_ind);

      // LIN_ALG::GEINV<dca::linalg::CPU>::execute_on_Green_function_matrix(f_matrix);
      geinv_obj.execute(f_matrix);

      for (int i = 0; i < N; i++)
        f_matrix(i, i) += I;

      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
          f_k_w(i, j, k_ind, w_ind) = f_matrix(i, j);
    }
  }
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_TRANSFORM_TO_ALPHA_HPP
