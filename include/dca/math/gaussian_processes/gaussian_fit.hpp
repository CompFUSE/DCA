// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a Gaussian fit.

#ifndef DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_FIT_HPP
#define DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_FIT_HPP

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"

namespace dca {
namespace math {
namespace gp {
// dca::math::gp::

template <typename scalar_type, typename K_dmn_t, typename k_dmn_t>
class gaussian_fit {
public:
  typedef typename K_dmn_t::parameter_type::element_type K_element_type;
  typedef typename k_dmn_t::parameter_type::element_type k_element_type;

  typedef linalg::Matrix<scalar_type, linalg::CPU> matrix_type;

  gaussian_fit();

  matrix_type& get_interpolation_matrix() {
    return interpolation_matrix;
  }
  matrix_type& get_correlation_matrix() {
    return correlation_matrix;
  }

  template <typename covariance_type>
  void initialize(covariance_type& covariance, double error);

private:
  template <typename covariance_type>
  void initialize_interpolation_matrix(covariance_type& covariance, double error);
  template <typename covariance_type>
  void initialize_correlation_matrix(covariance_type& covariance, double error);

  const static int DIMENSION = K_dmn_t::parameter_type::DIMENSION;

  std::vector<scalar_type> length_scale;

  matrix_type interpolation_matrix;
  matrix_type correlation_matrix;
};

template <typename scalar_type, typename K_dmn_t, typename k_dmn_t>
gaussian_fit<scalar_type, K_dmn_t, k_dmn_t>::gaussian_fit()
    : length_scale(DIMENSION, 1.),
      interpolation_matrix("interpolation_matrix"),
      correlation_matrix("correlation_matrix") {}

template <typename scalar_type, typename K_dmn_t, typename k_dmn_t>
template <typename covariance_type>
void gaussian_fit<scalar_type, K_dmn_t, k_dmn_t>::initialize(covariance_type& covariance,
                                                             double error) {
  initialize_interpolation_matrix(covariance, error);
  initialize_correlation_matrix(covariance, error);
}

template <typename scalar_type, typename K_dmn_t, typename k_dmn_t>
template <typename covariance_type>
void gaussian_fit<scalar_type, K_dmn_t, k_dmn_t>::initialize_interpolation_matrix(
    covariance_type& covariance, double error) {
  int N_r = k_dmn_t::dmn_size();
  int N_c = K_dmn_t::dmn_size();

  interpolation_matrix.resizeNoCopy(std::pair<int, int>(N_r, N_c));

  linalg::Matrix<scalar_type, linalg::CPU> K_K("K_K", std::pair<int, int>(N_c, N_c));

  for (int i = 0; i < N_c; i++) {
    for (int j = 0; j < N_c; j++) {
      K_element_type& x_i = K_dmn_t::get_elements()[i];
      K_element_type& x_j = K_dmn_t::get_elements()[j];

      std::vector<scalar_type> delta(DIMENSION, 0);
      for (int li = 0; li < DIMENSION; li++)
        delta[li] = (x_j[li] - x_i[li]);

      K_K(i, j) = covariance.execute(delta);  // compute_factor(periodic, delta);
    }

    K_K(i, i) += error * error;
  }

  linalg::Matrix<scalar_type, linalg::CPU> K_K_inv("K_K_inv", std::pair<int, int>(N_c, N_c));

  linalg::matrixop::pseudoInverse(K_K, K_K_inv);

  linalg::Matrix<scalar_type, linalg::CPU> interpolation_matrix_tmp("K_K_inv",
                                                                    std::pair<int, int>(N_r, N_c));

  for (int i = 0; i < N_r; i++) {
    for (int j = 0; j < N_c; j++) {
      k_element_type& x_i = k_dmn_t::get_elements()[i];
      K_element_type& x_j = K_dmn_t::get_elements()[j];

      std::vector<scalar_type> delta(DIMENSION, 0);
      for (int li = 0; li < DIMENSION; li++)
        delta[li] = (x_j[li] - x_i[li]);

      interpolation_matrix_tmp(i, j) = covariance.execute(delta);
    }

    interpolation_matrix_tmp(i, i) += error * error;
  }

  linalg::matrixop::gemm('N', 'N', interpolation_matrix_tmp, K_K_inv, interpolation_matrix);
}

template <typename scalar_type, typename K_dmn_t, typename k_dmn_t>
template <typename covariance_type>
void gaussian_fit<scalar_type, K_dmn_t, k_dmn_t>::initialize_correlation_matrix(
    covariance_type& covariance, double error) {
  int N_c = K_dmn_t::dmn_size();
  int N_r = k_dmn_t::dmn_size();

  correlation_matrix.resizeNoCopy(std::pair<int, int>(N_r, N_r));

  for (int i = 0; i < N_r; i++) {
    for (int j = 0; j < N_r; j++) {
      k_element_type& x_i = k_dmn_t::get_elements()[i];
      k_element_type& x_j = k_dmn_t::get_elements()[j];

      std::vector<scalar_type> delta(DIMENSION, 0);
      for (int li = 0; li < DIMENSION; li++)
        delta[li] = (x_j[li] - x_i[li]);

      correlation_matrix(i, j) = covariance.execute(delta);
    }

    correlation_matrix(i, i) += error * error;
  }

  linalg::Matrix<scalar_type, linalg::CPU> k_to_K("k_to_K", std::pair<int, int>(N_c, N_r));

  for (int i = 0; i < N_c; i++) {
    for (int j = 0; j < N_r; j++) {
      K_element_type& x_i = K_dmn_t::get_elements()[i];
      k_element_type& x_j = k_dmn_t::get_elements()[j];

      std::vector<scalar_type> delta(DIMENSION, 0);
      for (int li = 0; li < DIMENSION; li++)
        delta[li] = (x_j[li] - x_i[li]);

      correlation_matrix(i, j) = covariance.execute(delta);
    }

    k_to_K(i, i) += error * error;
  }

  linalg::matrixop::gemm('N', 'N', scalar_type(-1), interpolation_matrix, k_to_K, scalar_type(1),
                         correlation_matrix);
}

}  // gp
}  // math
}  // dca

#endif  // DCA_MATH_GAUSSIAN_PROCESSES_GAUSSIAN_FIT_HPP
