// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Rasmussen and Williams, p. 19, eq. 2.31.

#ifndef DCA_MATH_GAUSSIAN_PROCESSES_COVARIANCE_FUNCTION_HPP
#define DCA_MATH_GAUSSIAN_PROCESSES_COVARIANCE_FUNCTION_HPP

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "dca/linalg/matrix.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace math {
namespace gp {
// dca::math::gp::

enum GPCovariance { SQUARED_EXPONENTIAL, PERIODIC_SQUARED_EXPONENTIAL };

//
// Empty class template
//
template <GPCovariance covariance, typename k_dmn_t>
class covariance_function {};

//
// Partial specialization for squared exponential covariance
//
template <typename k_dmn_t>
class covariance_function<SQUARED_EXPONENTIAL, k_dmn_t> {
public:
  covariance_function();
  double execute(const std::vector<double>& x_i) const;
  void plot() const;

private:
  const static int DIMENSION = k_dmn_t::parameter_type::DIMENSION;

  double sigma_f;
  linalg::Matrix<double, linalg::CPU> A;
};

template <typename k_dmn_t>
covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::covariance_function()
    : sigma_f(1.), A("A", std::pair<int, int>(DIMENSION, DIMENSION)) {
  for (int i = 0; i < DIMENSION; i++)
    A(i, i) = 1.;
}

template <typename k_dmn_t>
double covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::execute(const std::vector<double>& x_i) const {
  std::vector<double> y_i(DIMENSION, 0.);

  for (int li = 0; li < DIMENSION; li++)
    for (int lj = 0; lj < DIMENSION; lj++)
      y_i[li] += k_dmn_t::parameter_type::get_inverse_basis()[li + lj * DIMENSION] * x_i[lj];

  double result = 0;

  for (int li = 0; li < DIMENSION; li++)
    for (int lj = 0; lj < DIMENSION; lj++)
      result += y_i[li] * A(li, lj) * y_i[lj];

  return sigma_f * std::exp(-result);
}

template <typename k_dmn_t>
void covariance_function<SQUARED_EXPONENTIAL, k_dmn_t>::plot() const {
  std::vector<double> x(0);
  std::vector<double> y(0);
  std::vector<double> z(0);

  std::vector<double> vec(DIMENSION, 0);

  switch (DIMENSION) {
    case 1: {
      for (int l = 0; l < k_dmn_t::dmn_size(); l++) {
        vec[0] = k_dmn_t::get_elements()[l];

        x.push_back(vec[0]);
        z.push_back(execute(vec));
      }

      util::Plot::plotPoints(x, z);
    } break;

    default:
      std::cout << __FUNCTION__ << std::endl;
  }
}

//
// Partial specialization for periodic squared exponential covariance
//
template <typename k_dmn_t>
class covariance_function<PERIODIC_SQUARED_EXPONENTIAL, k_dmn_t> {
public:
  covariance_function();
  double execute(const std::vector<double>& x_i) const;
  void plot() const;

private:
  const static int DIMENSION = k_dmn_t::parameter_type::DIMENSION;

  double sigma_f;
  linalg::Matrix<double, linalg::CPU> A;
};

template <typename k_dmn_t>
covariance_function<PERIODIC_SQUARED_EXPONENTIAL, k_dmn_t>::covariance_function()
    : sigma_f(1.), A("A", std::pair<int, int>(DIMENSION, DIMENSION)) {
  linalg::Matrix<double, linalg::CPU> T("T", std::pair<int, int>(DIMENSION, DIMENSION));

  for (int li = 0; li < DIMENSION; li++)
    for (int lj = 0; lj < DIMENSION; lj++)
      for (int lk = 0; lk < DIMENSION; lk++)
        T(li, lj) += k_dmn_t::parameter_type::get_inverse_super_basis()[li + lk * DIMENSION] *
                     k_dmn_t::parameter_type::get_basis()[lk + lj * DIMENSION];

  double l = 0;
  for (int lj = 0; lj < DIMENSION; lj++) {
    double norm = 0;
    for (int li = 0; li < DIMENSION; li++)
      norm += T(li, lj) * T(li, lj);

    l += sqrt(norm);
  }

  l /= DIMENSION;

  for (int li = 0; li < DIMENSION; li++)
    A(li, li) = 1. / (l * l);
}

template <typename k_dmn_t>
double covariance_function<PERIODIC_SQUARED_EXPONENTIAL, k_dmn_t>::execute(
    const std::vector<double>& x_i) const {
  std::vector<double> y_i(DIMENSION, 0.);

  for (int li = 0; li < DIMENSION; li++)
    for (int lj = 0; lj < DIMENSION; lj++)
      y_i[li] += k_dmn_t::parameter_type::get_inverse_super_basis()[li + lj * DIMENSION] * x_i[lj];

  for (int li = 0; li < DIMENSION; li++)
    y_i[li] = std::sin(M_PI * y_i[li]);

  double result = 0;

  for (int li = 0; li < DIMENSION; li++)
    for (int lj = 0; lj < DIMENSION; lj++)
      result += y_i[li] * A(li, lj) * y_i[lj];

  return sigma_f * std::exp(-result);
}

template <typename k_dmn_t>
void covariance_function<PERIODIC_SQUARED_EXPONENTIAL, k_dmn_t>::plot() const {
  double* basis = k_dmn_t::parameter_type::get_super_basis();

  std::vector<double> x(0);
  std::vector<double> y(0);
  std::vector<double> z(0);

  std::vector<double> vec(DIMENSION, 0);

  switch (DIMENSION) {
    case 1: {
      for (int i = -1; i <= 1; i++) {
        for (int l = 0; l < k_dmn_t::dmn_size(); l++) {
          vec[0] = k_dmn_t::get_elements()[l][0] + i * basis[0];

          x.push_back(vec[0]);

          z.push_back(execute(vec));
        }
      }

      util::Plot::plotPoints(x, z);
    } break;

    case 2: {
      for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
          for (int l = 0; l < k_dmn_t::dmn_size(); l++) {
            vec[0] = k_dmn_t::get_elements()[l][0] + i * basis[0] + j * basis[2];
            vec[1] = k_dmn_t::get_elements()[l][1] + i * basis[1] + j * basis[3];

            x.push_back(vec[0]);
            y.push_back(vec[1]);

            z.push_back(execute(vec));
          }
        }
      }

      util::Plot::heatMap(x, y, z);
    } break;

    default:
      std::cout << __FUNCTION__ << std::endl;
  }
}

}  // gp
}  // math
}  // dca

#endif  // DCA_MATH_GAUSSIAN_PROCESSES_COVARIANCE_FUNCTION_HPP
