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

namespace dca {
namespace math {
namespace inference {
// dca::math::inference::

template <typename DeconvolutionDmn, typename OtherDmn>
class RichardsonLucyDeconvolution {
public:
  RichardsonLucyDeconvolution(const double tolerance, const int max_iterations);

  // Returns the number of iterations executed.
  int findTargetFunction(
      const linalg::Matrix<double, linalg::CPU>& p,
      const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source,
      func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target);

  // Returns the number of iterations executed.
  int findTargetFunction(
      const linalg::Matrix<double, linalg::CPU>& p,
      const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source,
      func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target,
      func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target_convoluted);

private:
  void initializeMatrices(
      const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source);

  bool finished();

private:
  const double tolerance_;
  const int max_iterations_;

  func::function<bool, OtherDmn> is_finished_;

  linalg::Matrix<double, linalg::CPU> c_;
  linalg::Matrix<double, linalg::CPU> d_;
  linalg::Matrix<double, linalg::CPU> d_over_c_;
  linalg::Matrix<double, linalg::CPU> u_t_;
  linalg::Matrix<double, linalg::CPU> u_t_plus_1_;
};

template <typename DeconvolutionDmn, typename OtherDmn>
RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn>::RichardsonLucyDeconvolution(
    const double tolerance, const int max_iterations)
    : tolerance_(tolerance),
      max_iterations_(max_iterations),

      is_finished_("is_finished"),

      c_("c (Richardson-Lucy-deconvolution)"),
      d_("d (Richardson-Lucy-deconvolution)"),
      d_over_c_("d/c (Richardson-Lucy-deconvolution)"),
      u_t_("u_t (Richardson-Lucy-deconvolution)"),
      u_t_plus_1_("u_{t+1} (Richardson-Lucy-deconvolution)") {}

template <typename DeconvolutionDmn, typename OtherDmn>
int RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn>::findTargetFunction(
    const linalg::Matrix<double, linalg::CPU>& p,
    const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source,
    func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target) {
  assert(p.size().first == DeconvolutionDmn::dmn_size());
  assert(p.is_square());

  for (int i = 0; i < OtherDmn::dmn_size(); ++i)
    is_finished_(i) = false;

  initializeMatrices(source);

  int iterations = 0;
  while (!finished() && iterations < max_iterations_) {
    // Compute c.
    linalg::matrixop::gemm(p, u_t_, c_);

    // Compute d_over_c.
    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < DeconvolutionDmn::dmn_size(); ++i)
        d_over_c_(i, j) = d_(i, j) / c_(i, j);

    // Compute u_{t+1}.
    linalg::matrixop::gemm('T', 'N', p, d_over_c_, u_t_plus_1_);

    for (int j = 0; j < OtherDmn::dmn_size(); ++j)
      for (int i = 0; i < DeconvolutionDmn::dmn_size(); ++i)
        u_t_(i, j) = u_t_plus_1_(i, j) * u_t_(i, j);

    ++iterations;
  }

  // Copy iterative solution matrix into returned target function.
  for (int j = 0; j < OtherDmn::dmn_size(); ++j)
    for (int i = 0; i < DeconvolutionDmn::dmn_size(); ++i)
      target(i, j) = u_t_(i, j);

  return iterations;
}

template <typename DeconvolutionDmn, typename OtherDmn>
int RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn>::findTargetFunction(
    const linalg::Matrix<double, linalg::CPU>& p,
    const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source,
    func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target,
    func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& target_convoluted) {
  const int iterations = findTargetFunction(p, source, target);

  // Compute the convolution of the target function, which should resemble the source function.
  linalg::matrixop::gemm(p, u_t_, c_);

  for (int j = 0; j < OtherDmn::dmn_size(); j++)
    for (int i = 0; i < DeconvolutionDmn::dmn_size(); i++)
      target_convoluted(i, j) = c_(i, j);

  return iterations;
}

template <typename DeconvolutionDmn, typename OtherDmn>
void RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn>::initializeMatrices(
    const func::function<double, func::dmn_variadic<DeconvolutionDmn, OtherDmn>>& source) {
  const int num_rows = DeconvolutionDmn::dmn_size();
  const int num_cols = OtherDmn::dmn_size();

  c_.resizeNoCopy(std::make_pair(num_rows, num_cols));
  d_.resizeNoCopy(std::make_pair(num_rows, num_cols));
  d_over_c_.resizeNoCopy(std::make_pair(num_rows, num_cols));
  u_t_.resizeNoCopy(std::make_pair(num_rows, num_cols));
  u_t_plus_1_.resizeNoCopy(std::make_pair(num_rows, num_cols));

  // Initialize d matrix ("observed image") with source function.
  for (int j = 0; j < num_cols; ++j)
    for (int i = 0; i < num_rows; ++i)
      d_(i, j) = source(i, j);

  // Initialize iterative solution u_t with signs of column means of d.
  for (int j = 0; j < num_cols; ++j) {
    double mean = 0.;
    for (int i = 0; i < num_rows; ++i)
      mean += d_(i, j);
    mean /= num_rows;

    for (int i = 0; i < num_rows; ++i)
      u_t_(i, j) = mean / std::abs(mean);  // Used to be: u_t_(i,j) = mean.
  }

  // Initialize the other matrices with zero.
  for (int j = 0; j < num_cols; ++j)
    for (int i = 0; i < num_rows; ++i)
      c_(i, j) = d_over_c_(i, j) = u_t_plus_1_(i, j) = 0.;
}

template <typename DeconvolutionDmn, typename OtherDmn>
bool RichardsonLucyDeconvolution<DeconvolutionDmn, OtherDmn>::finished() {
  bool all_finished = true;

  for (int j = 0; j < OtherDmn::dmn_size(); ++j) {
    if (!is_finished_(j)) {
      double diff_squared = 0.;
      double norm_d_squared = 0.;

      // TODO: Fix error computation.
      for (int i = 0; i < DeconvolutionDmn::dmn_size(); ++i) {
        diff_squared += std::pow(c_(i, j) - d_(i, j), 2);
        norm_d_squared += std::pow(d_(i, j), 2);
      }

      const double error = std::sqrt(diff_squared / norm_d_squared);

      if (error < tolerance_)
        is_finished_(j) = true;

      else
        all_finished = false;
    }
  }

  return all_finished;
}

}  // inference
}  // math
}  // dca

#endif  // DCA_MATH_INFERENCE_RICHARDSON_LUCY_DECONVOLUTION_HPP
