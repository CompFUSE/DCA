// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the Richardson Lucy deconvolution algorithm.

#ifndef DCA_MATH_INFERENCE_RICHARDSON_LUCY_DECONVOLUTION_HPP
#define DCA_MATH_INFERENCE_RICHARDSON_LUCY_DECONVOLUTION_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace math {
namespace inference {
// dca::math::inference::

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
class RichardsonLucyDeconvolution {
public:
  RichardsonLucyDeconvolution(parameters_type& parameters_ref);

  void execute(const linalg::Matrix<double, linalg::CPU>& p,
               func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
               func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target);

  void execute(const linalg::Matrix<double, linalg::CPU>& p,
               func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
               func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_approx,
               func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target);

private:
  void initialize_matrices(func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source);

  void initialize_errors(func::function<bool, p_dmn_t>& is_finished,
                         func::function<double, p_dmn_t>& error_function);

  bool update_f_target(func::function<bool, p_dmn_t>& is_finished,
                       func::function<double, p_dmn_t>& error_function,
                       func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target);

private:
  parameters_type& parameters;
  typename parameters_type::concurrency_type& concurrency;

  linalg::Matrix<double, linalg::CPU> c;
  linalg::Matrix<double, linalg::CPU> d;

  linalg::Matrix<double, linalg::CPU> d_over_c;

  linalg::Matrix<double, linalg::CPU> u_t;
  linalg::Matrix<double, linalg::CPU> u_t_p_1;
};

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::RichardsonLucyDeconvolution(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      c("c (Richardson_Lucy_deconvolution)"),
      d("d (Richardson_Lucy_deconvolution)"),

      d_over_c("d/c (Richardson_Lucy_deconvolution)"),

      u_t("u_t (Richardson_Lucy_deconvolution)"),
      u_t_p_1("u_{t+1} (Richardson_Lucy_deconvolution)") {}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
void RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(
    const linalg::Matrix<double, linalg::CPU>& p,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target) {
  assert(p.size().first == k_dmn_t::dmn_size());
  assert(p.size().first == p.size().second);

  func::function<bool, p_dmn_t> is_finished("is_finished");
  func::function<double, p_dmn_t> error_function("error_function");

  initialize_matrices(f_source);

  // compute c
  linalg::matrixop::gemm(p, u_t, c);

  initialize_errors(is_finished, error_function);

  int l = 0;
  for (l = 0; l < parameters.get_deconvolution_iterations(); l++) {
    for (int j = 0; j < p_dmn_t::dmn_size(); j++)
      for (int i = 0; i < k_dmn_t::dmn_size(); i++)
        d_over_c(i, j) = d(i, j) / c(i, j);

    // compute u_t_plus_1
    linalg::matrixop::gemm('T', 'N', p, d_over_c, u_t_p_1);

    for (int j = 0; j < p_dmn_t::dmn_size(); j++)
      for (int i = 0; i < k_dmn_t::dmn_size(); i++)
        u_t(i, j) = u_t_p_1(i, j) * u_t(i, j);

    // compute c
    linalg::matrixop::gemm(p, u_t, c);

    bool finished = update_f_target(is_finished, error_function, f_target);

    if (finished)
      break;
  }

  for (int j = 0; j < p_dmn_t::dmn_size(); j++)
    if (not is_finished(j))
      for (int i = 0; i < k_dmn_t::dmn_size(); i++)
        f_target(i, j) = u_t(i, j);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t\t Richardson-Lucy deconvolution: " << l << " iterations" << std::endl;
  }
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
void RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(
    const linalg::Matrix<double, linalg::CPU>& p,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_approx,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target) {
  execute(p, f_source, f_target);

  for (int j = 0; j < p_dmn_t::dmn_size(); j++)
    for (int i = 0; i < k_dmn_t::dmn_size(); i++)
      u_t(i, j) = f_target(i, j);

  linalg::matrixop::gemm(p, u_t, c);

  for (int j = 0; j < p_dmn_t::dmn_size(); j++)
    for (int i = 0; i < k_dmn_t::dmn_size(); i++)
      f_approx(i, j) = c(i, j);
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
void RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_matrices(
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source) {
  int nr_rows = k_dmn_t::dmn_size();
  int nr_cols = p_dmn_t::dmn_size();

  c.resizeNoCopy(std::pair<int, int>(nr_rows, nr_cols));
  d.resizeNoCopy(std::pair<int, int>(nr_rows, nr_cols));

  d_over_c.resizeNoCopy(std::pair<int, int>(nr_rows, nr_cols));

  u_t.resizeNoCopy(std::pair<int, int>(nr_rows, nr_cols));
  u_t_p_1.resizeNoCopy(std::pair<int, int>(nr_rows, nr_cols));

  for (int j = 0; j < nr_cols; j++)
    for (int i = 0; i < nr_rows; i++)
      d(i, j) = f_source(i, j);

  for (int j = 0; j < nr_cols; j++) {
    double mean = 0;
    for (int i = 0; i < nr_rows; i++)
      mean += f_source(i, j);
    mean /= nr_rows;

    for (int i = 0; i < nr_rows; i++)
      u_t(i, j) = mean / std::abs(mean);  // u_t(i,j) = mean;
  }
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
void RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::initialize_errors(
    func::function<bool, p_dmn_t>& is_finished, func::function<double, p_dmn_t>& error_function) {
  for (int j = 0; j < p_dmn_t::dmn_size(); j++) {
    double error = 0;
    for (int i = 0; i < k_dmn_t::dmn_size(); i++)
      error += std::pow(c(i, j) - d(i, j), 2);

    is_finished(j) = false;
    error_function(j) = sqrt(error) / k_dmn_t::dmn_size();
  }
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
bool RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::update_f_target(
    func::function<bool, p_dmn_t>& is_finished, func::function<double, p_dmn_t>& error_function,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target) {
  bool all_are_finished = true;

  double epsilon = parameters.get_deconvolution_tolerance();

  for (int j = 0; j < p_dmn_t::dmn_size(); j++) {
    if (not is_finished(j)) {
      double diff = 0;
      double tot = 1.e-6;

      for (int i = 0; i < k_dmn_t::dmn_size(); i++) {
        diff += std::pow(c(i, j) - d(i, j), 2);
        tot += std::pow(d(i, j), 2);
      }

      error_function(j) = std::sqrt(diff / tot);

      if (error_function(j) < epsilon) {
        for (int i = 0; i < k_dmn_t::dmn_size(); i++)
          f_target(i, j) = u_t(i, j);

        is_finished(j) = true;
      }

      else {
        all_are_finished = false;
      }
    }
  }

  return all_are_finished;
}

}  // inference
}  // math
}  // dca

#endif  // DCA_MATH_INFERENCE_RICHARDSON_LUCY_DECONVOLUTION_HPP
