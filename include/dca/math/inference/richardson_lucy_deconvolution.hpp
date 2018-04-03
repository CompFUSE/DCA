// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
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

  // Returns the number of iterations executed.
  int execute(const linalg::Matrix<double, linalg::CPU>& p,
              const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
              func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target);

  // Returns the number of iterations executed.
  int execute(const linalg::Matrix<double, linalg::CPU>& p,
              const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
              func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_approx,
              func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target);

private:
  void initializeMatrices(const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& source);

  bool computeError(func::function<bool, p_dmn_t>& is_finished);

private:
  parameters_type& parameters;

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

      c("c (Richardson_Lucy_deconvolution)"),
      d("d (Richardson_Lucy_deconvolution)"),

      d_over_c("d/c (Richardson_Lucy_deconvolution)"),

      u_t("u_t (Richardson_Lucy_deconvolution)"),
      u_t_p_1("u_{t+1} (Richardson_Lucy_deconvolution)") {}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
int RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(
    const linalg::Matrix<double, linalg::CPU>& p,
    const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target) {
  assert(p.size().first == k_dmn_t::dmn_size());
  assert(p.size().first == p.size().second);

  func::function<bool, p_dmn_t> is_finished("is_finished");
  for (int i = 0; i < p_dmn_t::dmn_size(); ++i)
    is_finished(i) = false;

  initializeMatrices(f_source);

  // compute c
  linalg::matrixop::gemm(p, u_t, c);

  int iters = 0;
  for (; iters < parameters.get_deconvolution_iterations(); ++iters) {
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

    bool finished = computeError(is_finished);

    if (finished)
      break;
  }

  // Copy iterative solution matrix into returned target function.
  for (int j = 0; j < p_dmn_t::dmn_size(); ++j)
    for (int i = 0; i < k_dmn_t::dmn_size(); ++i)
      f_target(i, j) = u_t(i, j);

  return iters;
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
int RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::execute(
    const linalg::Matrix<double, linalg::CPU>& p,
    const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_source,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_approx,
    func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& f_target) {
  const int iterations = execute(p, f_source, f_target);

  for (int j = 0; j < p_dmn_t::dmn_size(); j++)
    for (int i = 0; i < k_dmn_t::dmn_size(); i++)
      u_t(i, j) = f_target(i, j);

  linalg::matrixop::gemm(p, u_t, c);

  for (int j = 0; j < p_dmn_t::dmn_size(); j++)
    for (int i = 0; i < k_dmn_t::dmn_size(); i++)
      f_approx(i, j) = c(i, j);

  return iterations;
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
void RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::initializeMatrices(
    const func::function<double, func::dmn_variadic<k_dmn_t, p_dmn_t>>& source) {
  const int num_rows = k_dmn_t::dmn_size();
  const int num_cols = p_dmn_t::dmn_size();

  c.resizeNoCopy(std::make_pair(num_rows, num_cols));
  d.resizeNoCopy(std::make_pair(num_rows, num_cols));
  d_over_c.resizeNoCopy(std::make_pair(num_rows, num_cols));
  u_t.resizeNoCopy(std::make_pair(num_rows, num_cols));
  u_t_p_1.resizeNoCopy(std::make_pair(num_rows, num_cols));

  // Initialize "observed image" d (source).
  for (int j = 0; j < num_cols; ++j)
    for (int i = 0; i < num_rows; ++i)
      d(i, j) = source(i, j);

  // Initialize iterative solution u_t (target) with sign of column means of d.
  for (int j = 0; j < num_cols; ++j) {
    double mean = 0.;
    for (int i = 0; i < num_rows; ++i)
      mean += d(i, j);
    mean /= num_rows;

    for (int i = 0; i < num_rows; ++i)
      u_t(i, j) = mean / std::abs(mean);  // Used to be: u_t(i,j) = mean.
  }
}

template <typename parameters_type, typename k_dmn_t, typename p_dmn_t>
bool RichardsonLucyDeconvolution<parameters_type, k_dmn_t, p_dmn_t>::computeError(
    func::function<bool, p_dmn_t>& is_finished) {
  bool all_are_finished = true;

  const double epsilon = parameters.get_deconvolution_tolerance();

  for (int j = 0; j < p_dmn_t::dmn_size(); ++j) {
    if (!is_finished(j)) {
      double diff_squared = 0.;
      double norm_d_squared = 0.;

      // TODO: Fix error computation.
      for (int i = 0; i < k_dmn_t::dmn_size(); ++i) {
        diff_squared += std::pow(c(i, j) - d(i, j), 2);
        norm_d_squared += std::pow(d(i, j), 2);
      }

      const double error = std::sqrt(diff_squared / norm_d_squared);

      if (error < epsilon)
        is_finished(j) = true;

      else
        all_are_finished = false;
    }
  }

  return all_are_finished;
}

}  // inference
}  // math
}  // dca

#endif  // DCA_MATH_INFERENCE_RICHARDSON_LUCY_DECONVOLUTION_HPP
