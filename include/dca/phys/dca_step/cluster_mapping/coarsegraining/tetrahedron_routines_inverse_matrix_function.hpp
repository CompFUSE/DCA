// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Tetrahedron routines: inverse matrix function.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_INVERSE_MATRIX_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_INVERSE_MATRIX_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_eigenvalue_degeneracy.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_integration_data.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

class tetrahedron_routines_inverse_matrix_function {
  inline static double EPSILON() {
    return 1.e-2;
  }

  template <typename scalartype>
  struct matrix_element_struct {
    int i;

    std::complex<scalartype> e;
    std::complex<scalartype> f;

    std::complex<scalartype> log_min_e;

    std::complex<scalartype> r;
  };

private:
  // TODO: Replace this method with a more efficient one.
  template <typename scalartype>
  inline static std::complex<scalartype> Power(std::complex<scalartype> val, int n);

  template <typename scalartype>
  static bool are_equal(std::complex<scalartype> const& x, std::complex<scalartype> const& y);

  template <typename scalartype>
  static bool not_equal(std::complex<scalartype> const& x, std::complex<scalartype> const& y);

  template <typename scalartype>
  static bool value_comp(matrix_element_struct<scalartype> const& x,
                         matrix_element_struct<scalartype> const& y);

  template <typename scalartype>
  static bool permutation_comp(matrix_element_struct<scalartype> const& x,
                               matrix_element_struct<scalartype> const& y);

public:
  // 1D
  template <typename scalartype>
  static void execute(int size, scalartype volume, std::complex<scalartype>* G_0,
                      std::complex<scalartype>* G_1, std::complex<scalartype>* f_result,
                      tetrahedron_integration_data<scalartype>& data_obj);

  // 2D
  template <typename scalartype>
  static void execute(int size, scalartype volume, std::complex<scalartype>* G_0,
                      std::complex<scalartype>* G_1, std::complex<scalartype>* G_2,
                      std::complex<scalartype>* f_result,
                      tetrahedron_integration_data<scalartype>& data_obj);

  // 3D
  template <typename scalartype>
  static void execute(int size, scalartype volume, std::complex<scalartype>* G_0,
                      std::complex<scalartype>* G_1, std::complex<scalartype>* G_2,
                      std::complex<scalartype>* G_3, std::complex<scalartype>* f_result,
                      tetrahedron_integration_data<scalartype>& data_obj);

private:
  template <typename scalartype>
  static std::complex<scalartype> integrate_matrix_element_1D(
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static std::complex<scalartype> integrate_matrix_element_2D(
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static void integrate_eigenvalues_2D(math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
                                       std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static std::complex<scalartype> integrate_matrix_element_3D(
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static std::complex<scalartype> integrate_matrix_element_3D(
      math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static void integrate_eigenvalues_3D(math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
                                       std::vector<matrix_element_struct<scalartype>>& vec);
  template <typename scalartype>
  static math::geometry::TetrahedronEigenvalueDegeneracy find_degeneracy_1D(
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static math::geometry::TetrahedronEigenvalueDegeneracy find_degeneracy_2D(
      std::vector<matrix_element_struct<scalartype>>& vec);

  template <typename scalartype>
  static math::geometry::TetrahedronEigenvalueDegeneracy find_degeneracy_3D(
      std::vector<matrix_element_struct<scalartype>>& vec);
};

template <typename scalartype>
std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::Power(
    std::complex<scalartype> val, int n) {
  switch (n) {
    case 1:
      return val;
      break;

    case 2:
      return val * val;
      break;

    case 3:
      return val * val * val;
      break;

    case 4:
      return val * val * val * val;
      break;

    case 5:
      return val * val * val * val * val;
      break;

    default:
      return std::pow(val, n);
  }
}

template <typename scalartype>
bool tetrahedron_routines_inverse_matrix_function::are_equal(std::complex<scalartype> const& x,
                                                             std::complex<scalartype> const& y) {
  scalartype norm_x = std::norm(x);
  scalartype norm_y = std::norm(y);

  scalartype max = norm_x > norm_y ? norm_x : norm_y;

  max += 1.e-16;

  scalartype abs_error = std::norm(x - y);
  scalartype rel_error = abs_error / max;

  if (abs_error < EPSILON() * EPSILON() or rel_error < EPSILON() * EPSILON())
    return true;

  return false;
}

template <typename scalartype>
bool tetrahedron_routines_inverse_matrix_function::not_equal(std::complex<scalartype> const& x,
                                                             std::complex<scalartype> const& y) {
  return (not are_equal(x, y));
}
template <typename scalartype>
bool tetrahedron_routines_inverse_matrix_function::value_comp(
    matrix_element_struct<scalartype> const& x, matrix_element_struct<scalartype> const& y) {
  return std::norm(x.e) > std::norm(y.e);
}

template <typename scalartype>
bool tetrahedron_routines_inverse_matrix_function::permutation_comp(
    matrix_element_struct<scalartype> const& x, matrix_element_struct<scalartype> const& y) {
  return x.i < y.i;
}

/************************************
 ***
 ***   1D tetrahedron-integration
 ***
 ************************************/

template <typename scalartype>
void tetrahedron_routines_inverse_matrix_function::execute(
    int N, scalartype volume, std::complex<scalartype>* G_0, std::complex<scalartype>* G_1,
    std::complex<scalartype>* f_result, tetrahedron_integration_data<scalartype>& data_obj) {
  // diagonolize the G-matrices

  memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.G_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N,
                            data_obj.VR_0, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N,
                            data_obj.VR_1, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);

  memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.VR_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  // integrate G-matrices

  for (int l = 0; l < N * N; l++)
    f_result[l] = 0;

  std::vector<matrix_element_struct<scalartype>> vec(2);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int l = 0; l < N; l++) {
        vec[0].i = 0;
        vec[1].i = 1;

        vec[0].e = data_obj.W_0[l];
        vec[1].e = data_obj.W_1[l];

        vec[0].f = data_obj.VR_0[i + l * N] * data_obj.VR_inv_0[l + j * N];
        vec[1].f = data_obj.VR_1[i + l * N] * data_obj.VR_inv_1[l + j * N];

        f_result[i + j * N] += integrate_matrix_element_1D(vec);
      }
    }
  }

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      f_result[i + j * N] *= volume;
}

template <typename scalartype>
std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_1D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 2);

  math::geometry::TetrahedronEigenvalueDegeneracy degeneracy = find_degeneracy_1D(vec);

  std::complex<scalartype> r[2];
  std::complex<scalartype> f[3];

  f[0] = vec[0].f;
  f[1] = vec[1].f;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      r[0] = (e0 - e1 - e1 * std::log(-e0) + e1 * std::log(-e1)) / Power(e0 - e1, 2);
      r[1] = (-e0 + e1 + e0 * std::log(-e0) - e0 * std::log(-e1)) / Power(e0 - e1, 2);
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      r[0] = 1 / (2. * e0);
      r[1] = 1 / (2. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  std::complex<scalartype> result = 0;

  {
    result += f[0] * r[0];
    result += f[1] * r[1];

    assert(result == result);  // make sure there is no NAN !
  }

  return result;
}

template <typename scalartype>
math::geometry::TetrahedronEigenvalueDegeneracy tetrahedron_routines_inverse_matrix_function::find_degeneracy_1D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 2);

  if (are_equal(vec[0].e, vec[1].e))
    return math::geometry::TWOFOLD_DEGENERACY;

  return math::geometry::NO_DEGENERACY;
}

/************************************
 ***
 ***   2D tetrahedron-integration
 ***
 ************************************/

template <typename scalartype>
void tetrahedron_routines_inverse_matrix_function::execute(
    int N, scalartype volume, std::complex<scalartype>* G_0, std::complex<scalartype>* G_1,
    std::complex<scalartype>* G_2, std::complex<scalartype>* f_result,
    tetrahedron_integration_data<scalartype>& data_obj) {
  // diagonolize the G-matrices
  memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.G_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_2, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N,
                            data_obj.VR_0, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N,
                            data_obj.VR_1, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N,
                            data_obj.VR_2, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);

  memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.VR_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_2, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  for (int l = 0; l < N * N; l++)
    f_result[l] = 0;

  std::vector<matrix_element_struct<scalartype>> vec(3);

  if (false) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int l = 0; l < N; l++) {
          {
            vec[0].i = 0;
            vec[1].i = 1;
            vec[2].i = 2;

            vec[0].e = data_obj.W_0[l];
            vec[1].e = data_obj.W_1[l];
            vec[2].e = data_obj.W_2[l];

            vec[0].f = data_obj.VR_0[i + l * N] * data_obj.VR_inv_0[l + j * N];
            vec[1].f = data_obj.VR_1[i + l * N] * data_obj.VR_inv_1[l + j * N];
            vec[2].f = data_obj.VR_2[i + l * N] * data_obj.VR_inv_2[l + j * N];

            f_result[i + j * N] += integrate_matrix_element_2D(vec);
          }
        }
      }
    }
  }
  else {
    for (int l = 0; l < N; l++) {
      vec[0].i = 0;
      vec[1].i = 1;
      vec[2].i = 2;

      vec[0].e = data_obj.W_0[l];
      vec[1].e = data_obj.W_1[l];
      vec[2].e = data_obj.W_2[l];

      vec[0].log_min_e = std::log(-vec[0].e);
      vec[1].log_min_e = std::log(-vec[1].e);
      vec[2].log_min_e = std::log(-vec[2].e);

      math::geometry::TetrahedronEigenvalueDegeneracy degeneracy = find_degeneracy_2D(vec);

      integrate_eigenvalues_2D(degeneracy, vec);

      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          vec[0].f = data_obj.VR_0[i + l * N] * data_obj.VR_inv_0[l + j * N];
          vec[1].f = data_obj.VR_1[i + l * N] * data_obj.VR_inv_1[l + j * N];
          vec[2].f = data_obj.VR_2[i + l * N] * data_obj.VR_inv_2[l + j * N];

          f_result[i + j * N] += vec[0].r * vec[0].f;
          f_result[i + j * N] += vec[1].r * vec[1].f;
          f_result[i + j * N] += vec[2].r * vec[2].f;
        }
      }
    }
  }

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      f_result[i + j * N] *= (2. * volume);
}

template <typename scalartype>
std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_2D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 3);

  math::geometry::TetrahedronEigenvalueDegeneracy degeneracy = find_degeneracy_2D(vec);

  std::complex<scalartype> r[3];
  std::complex<scalartype> f[3];

  f[0] = vec[0].f;
  f[1] = vec[1].f;
  f[2] = vec[2].f;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;
  std::complex<scalartype> e2 = vec[2].e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      assert(not_equal(vec[0].e, vec[1].e));
      assert(not_equal(vec[1].e, vec[2].e));
      assert(not_equal(vec[2].e, vec[0].e));

      r[0] = (-(e0 * (e1 - e2) * (-2. * e1 * e2 + e0 * (e1 + e2)) * std::log(-e0)) +
              Power(e1, 2) * Power(e0 - e2, 2) * std::log(-e1) +
              (e0 - e1) * (e0 * (e0 - e2) * (e1 - e2) + (-e0 + e1) * Power(e2, 2) * std::log(-e2))) /
             (2. * Power(e0 - e1, 2) * Power(e0 - e2, 2) * (e1 - e2));
      r[1] = (Power(e0, 2) * Power(e1 - e2, 2) * std::log(-e0) +
              e1 * (e0 - e2) *
                  (-((e0 - e1) * (e1 - e2)) - (e0 * e1 - 2. * e0 * e2 + e1 * e2) * std::log(-e1)) -
              Power(e0 - e1, 2) * Power(e2, 2) * std::log(-e2)) /
             (2. * Power(e0 - e1, 2) * (e0 - e2) * Power(e1 - e2, 2));
      r[2] = (Power(e0, 2) * Power(e1 - e2, 2) * std::log(-e0) -
              Power(e1, 2) * Power(e0 - e2, 2) * std::log(-e1) +
              (-e0 + e1) * e2 *
                  ((e0 - e2) * (-e1 + e2) + (-2. * e0 * e1 + (e0 + e1) * e2) * std::log(-e2))) /
             (2. * (e0 - e1) * Power(e0 - e2, 2) * Power(e1 - e2, 2));
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      assert(not_equal(vec[0].e, vec[1].e));
      assert(are_equal(vec[1].e, vec[2].e));

      r[0] = (Power(e0, 2) - Power(e1, 2) - 2. * e0 * e1 * std::log(-e0) +
              2. * e0 * e1 * std::log(-e1)) /
             (2. * Power(e0 - e1, 3));
      r[1] = -(3. * Power(e0, 2) - 4. * e0 * e1 + Power(e1, 2) +
               2. * Power(e0, 2) * (-std::log(-e0) + std::log(-e1))) /
             (4. * Power(e0 - e1, 3));
      r[2] = -(3. * Power(e0, 2) - 4. * e0 * e1 + Power(e1, 2) +
               2. * Power(e0, 2) * (-std::log(-e0) + std::log(-e1))) /
             (4. * Power(e0 - e1, 3));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY: {
      assert(are_equal(vec[0].e, vec[1].e));
      assert(are_equal(vec[0].e, vec[2].e));

      r[0] = 1. / (6. * e0);
      r[1] = 1. / (6. * e0);
      r[2] = 1. / (6. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  std::complex<scalartype> result = 0;

  {
    result += f[0] * r[0];
    result += f[1] * r[1];
    result += f[2] * r[2];

    assert(result == result);  // make sure there is no NAN !
  }

  return result;
}

template <typename scalartype>
void tetrahedron_routines_inverse_matrix_function::integrate_eigenvalues_2D(
    math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 3);

  int i0 = vec[0].i;
  int i1 = vec[1].i;
  int i2 = vec[2].i;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;
  std::complex<scalartype> e2 = vec[2].e;

  std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
  std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
  std::complex<scalartype> log_min_e2 = vec[2].log_min_e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      assert(not_equal(vec[0].e, vec[1].e));
      assert(not_equal(vec[1].e, vec[2].e));
      assert(not_equal(vec[2].e, vec[0].e));

      vec[i0].r =
          (-(e0 * (e1 - e2) * (-2. * e1 * e2 + e0 * (e1 + e2)) * log_min_e0) +
           Power(e1, 2) * Power(e0 - e2, 2) * log_min_e1 +
           (e0 - e1) * (e0 * (e0 - e2) * (e1 - e2) + (-e0 + e1) * Power(e2, 2) * log_min_e2)) /
          (2. * Power(e0 - e1, 2) * Power(e0 - e2, 2) * (e1 - e2));
      vec[i1].r = (Power(e0, 2) * Power(e1 - e2, 2) * log_min_e0 +
                   e1 * (e0 - e2) *
                       (-((e0 - e1) * (e1 - e2)) - (e0 * e1 - 2. * e0 * e2 + e1 * e2) * log_min_e1) -
                   Power(e0 - e1, 2) * Power(e2, 2) * log_min_e2) /
                  (2. * Power(e0 - e1, 2) * (e0 - e2) * Power(e1 - e2, 2));
      vec[i2].r = (Power(e0, 2) * Power(e1 - e2, 2) * log_min_e0 -
                   Power(e1, 2) * Power(e0 - e2, 2) * log_min_e1 +
                   (-e0 + e1) * e2 *
                       ((e0 - e2) * (-e1 + e2) + (-2. * e0 * e1 + (e0 + e1) * e2) * log_min_e2)) /
                  (2. * (e0 - e1) * Power(e0 - e2, 2) * Power(e1 - e2, 2));
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      assert(not_equal(vec[0].e, vec[1].e));
      assert(are_equal(vec[1].e, vec[2].e));

      vec[i0].r =
          (Power(e0, 2) - Power(e1, 2) - 2. * e0 * e1 * log_min_e0 + 2. * e0 * e1 * log_min_e1) /
          (2. * Power(e0 - e1, 3));
      vec[i1].r = -(3. * Power(e0, 2) - 4. * e0 * e1 + Power(e1, 2) +
                    2. * Power(e0, 2) * (-log_min_e0 + log_min_e1)) /
                  (4. * Power(e0 - e1, 3));
      vec[i2].r = -(3. * Power(e0, 2) - 4. * e0 * e1 + Power(e1, 2) +
                    2. * Power(e0, 2) * (-log_min_e0 + log_min_e1)) /
                  (4. * Power(e0 - e1, 3));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY: {
      assert(are_equal(vec[0].e, vec[1].e));
      assert(are_equal(vec[0].e, vec[2].e));

      vec[i0].r = 1. / (6. * e0);
      vec[i1].r = 1. / (6. * e0);
      vec[i2].r = 1. / (6. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalartype>
math::geometry::TetrahedronEigenvalueDegeneracy tetrahedron_routines_inverse_matrix_function::find_degeneracy_2D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 3);

  if (not_equal(vec[0].e, vec[1].e) and not_equal(vec[1].e, vec[2].e) and
      not_equal(vec[2].e, vec[0].e))
    return math::geometry::NO_DEGENERACY;

  if (are_equal(vec[0].e, vec[1].e) and are_equal(vec[0].e, vec[2].e))
    return math::geometry::THREEFOLD_DEGENERACY;

  do {
    if (not_equal(vec[0].e, vec[1].e) and are_equal(vec[1].e, vec[2].e))
      return math::geometry::TWOFOLD_DEGENERACY;
  } while (std::next_permutation(
      vec.begin(), vec.end(),
      tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>));

  throw std::logic_error(__FUNCTION__);
  return math::geometry::NO_DEGENERACY;
}

/************************************
 ***
 ***   3D tetrahedron-integration
 ***
 ************************************/

template <typename scalartype>
void tetrahedron_routines_inverse_matrix_function::execute(
    int N, scalartype volume, std::complex<scalartype>* G_0, std::complex<scalartype>* G_1,
    std::complex<scalartype>* G_2, std::complex<scalartype>* G_3,
    std::complex<scalartype>* f_result, tetrahedron_integration_data<scalartype>& data_obj) {
  // obtain G^{-1} from G
  memcpy(data_obj.G_inv_0, G_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_1, G_1, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_2, G_2, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.G_inv_3, G_3, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.G_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_2, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.G_inv_3, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  // diagonolize the G-matrices
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_0, N, data_obj.W_0, data_obj.VR_inv_0, N,
                            data_obj.VR_0, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_1, N, data_obj.W_1, data_obj.VR_inv_1, N,
                            data_obj.VR_1, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_2, N, data_obj.W_2, data_obj.VR_inv_2, N,
                            data_obj.VR_2, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);
  dca::linalg::lapack::geev("N", "V", N, data_obj.G_inv_3, N, data_obj.W_3, data_obj.VR_inv_3, N,
                            data_obj.VR_3, N, data_obj.eig_work, data_obj.eig_lwork,
                            data_obj.eig_rwork);

  // obtain V^{-1}
  memcpy(data_obj.VR_inv_0, data_obj.VR_0, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_1, data_obj.VR_1, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_2, data_obj.VR_2, sizeof(std::complex<scalartype>) * N * N);
  memcpy(data_obj.VR_inv_3, data_obj.VR_3, sizeof(std::complex<scalartype>) * N * N);

  dca::linalg::lapack::inverse(N, data_obj.VR_inv_0, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_1, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_2, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);
  dca::linalg::lapack::inverse(N, data_obj.VR_inv_3, data_obj.inv_ipiv, data_obj.inv_work,
                               data_obj.inv_lwork);

  // integrate G-matrices
  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      f_result[i + j * N] = 0.;

  std::vector<matrix_element_struct<scalartype>> vec(4);

  for (int l = 0; l < N; l++) {
    vec[0].i = 0;
    vec[1].i = 1;
    vec[2].i = 2;
    vec[3].i = 3;

    vec[0].e = data_obj.W_0[l];
    vec[1].e = data_obj.W_1[l];
    vec[2].e = data_obj.W_2[l];
    vec[3].e = data_obj.W_3[l];

    vec[0].log_min_e = std::log(-vec[0].e);
    vec[1].log_min_e = std::log(-vec[1].e);
    vec[2].log_min_e = std::log(-vec[2].e);
    vec[3].log_min_e = std::log(-vec[3].e);

    math::geometry::TetrahedronEigenvalueDegeneracy degeneracy = find_degeneracy_3D(vec);

    integrate_eigenvalues_3D(degeneracy, vec);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        vec[0].f = data_obj.VR_0[i + l * N] * data_obj.VR_inv_0[l + j * N];
        vec[1].f = data_obj.VR_1[i + l * N] * data_obj.VR_inv_1[l + j * N];
        vec[2].f = data_obj.VR_2[i + l * N] * data_obj.VR_inv_2[l + j * N];
        vec[3].f = data_obj.VR_3[i + l * N] * data_obj.VR_inv_3[l + j * N];

        f_result[i + j * N] += vec[0].r * vec[0].f;
        f_result[i + j * N] += vec[1].r * vec[1].f;
        f_result[i + j * N] += vec[2].r * vec[2].f;
        f_result[i + j * N] += vec[3].r * vec[3].f;
      }
    }
  }

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      f_result[i + j * N] *= (6. * volume);
}

template <typename scalartype>
std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 4);

  math::geometry::TetrahedronEigenvalueDegeneracy degeneracy = find_degeneracy_3D(vec);

  std::complex<scalartype> r[4];
  std::complex<scalartype> f[4];

  f[0] = vec[0].f;
  f[1] = vec[1].f;
  f[2] = vec[2].f;
  f[3] = vec[3].f;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;
  std::complex<scalartype> e2 = vec[2].e;
  std::complex<scalartype> e3 = vec[3].e;

  std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
  std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
  std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
  std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e3));
      assert(not_equal(e3, e0));
      assert(not_equal(e0, e2));
      assert(not_equal(e1, e3));

      r[0] = (-((Power(e0, 2) * (3. * e1 * e2 * e3 + Power(e0, 2) * (e1 + e2 + e3) -
                                 2. * e0 * (e1 * e2 + (e1 + e2) * e3)) *
                 log_min_e0) /
                (Power(e0 - e1, 2) * Power(e0 - e2, 2) * Power(e0 - e3, 2))) +
              (Power(e1, 3) * log_min_e1) / (Power(e0 - e1, 2) * (e1 - e2) * (e1 - e3)) +
              (Power(e2, 3) * log_min_e2) / (Power(e0 - e2, 2) * (-e1 + e2) * (e2 - e3)) +
              ((Power(e0, 2) * (e0 - e3)) / ((e0 - e1) * (e0 - e2)) +
               (Power(e3, 3) * log_min_e3) / ((-e1 + e3) * (-e2 + e3))) /
                  Power(e0 - e3, 2)) /
             6.;
      r[1] = (Power(e1, 2) / ((-e0 + e1) * (e1 - e2) * (e1 - e3)) +
              (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * (e0 - e2) * (e0 - e3)) -
              (Power(e1, 2) * (e0 * (Power(e1, 2) + 3. * e2 * e3 - 2. * e1 * (e2 + e3)) +
                               e1 * (-2. * e2 * e3 + e1 * (e2 + e3))) *
               log_min_e1) /
                  (Power(e0 - e1, 2) * Power(e1 - e2, 2) * Power(e1 - e3, 2)) +
              ((Power(e2, 3) * log_min_e2) / (Power(e1 - e2, 2) * (-e0 + e2)) +
               (Power(e3, 3) * log_min_e3) / ((e0 - e3) * Power(e1 - e3, 2))) /
                  (e2 - e3)) /
             6.;
      r[2] = (Power(e0, 3) * log_min_e0 +
              (-(Power(e1, 3) * Power(e0 - e2, 2) * (e0 - e3) * Power(e2 - e3, 2) * log_min_e1) +
               (e0 - e1) *
                   (Power(e2, 2) * (e0 - e3) * (-e1 + e3) *
                        (e0 * (-2. * e1 * e2 + Power(e2, 2) + 3. * e1 * e3 - 2. * e2 * e3) +
                         e2 * (e1 * e2 - 2. * e1 * e3 + e2 * e3)) *
                        log_min_e2 +
                    (-e0 + e2) * (-e1 + e2) * (Power(e2, 2) * (e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                               (e0 - e2) * (e1 - e2) * Power(e3, 3) * log_min_e3))) /
                  (Power(e1 - e2, 2) * (e1 - e3) * Power(e2 - e3, 2))) /
             (6. * (e0 - e1) * Power(e0 - e2, 2) * (e0 - e3));
      r[3] = (Power(e0, 3) * log_min_e0 +
              (-(Power(e1, 3) * (e0 - e2) * Power(e0 - e3, 2) * Power(e2 - e3, 2) * log_min_e1) +
               (e0 - e1) * (Power(e2, 3) * Power(e0 - e3, 2) * Power(e1 - e3, 2) * log_min_e2 +
                            (e0 - e2) * (-e1 + e2) * Power(e3, 2) *
                                ((e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                 (3. * e0 * e1 * e2 - 2. * (e0 * e1 + (e0 + e1) * e2) * e3 +
                                  (e0 + e1 + e2) * Power(e3, 2)) *
                                     log_min_e3))) /
                  ((e1 - e2) * Power(e1 - e3, 2) * Power(e2 - e3, 2))) /
             (6. * (e0 - e1) * (e0 - e2) * Power(e0 - e3, 2));
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e0));
      assert(are_equal(e2, e3));

      r[0] = ((Power(e0, 2) * (e1 - e2) - e0 * Power(e2, 2) + e1 * Power(e2, 2)) /
                  ((e0 - e1) * Power(e0 - e2, 2) * (e1 - e2)) -
              (Power(e0, 2) * (-3. * e1 * e2 + e0 * (e1 + 2. * e2)) * log_min_e0) /
                  (Power(e0 - e1, 2) * Power(e0 - e2, 3)) +
              ((Power(e1, 3) * log_min_e1) / Power(e0 - e1, 2) +
               (Power(e2, 2) * (-3. * e0 * e1 + 2. * e0 * e2 + e1 * e2) * log_min_e2) /
                   Power(e0 - e2, 3)) /
                  Power(e1 - e2, 2)) /
             6.;
      r[1] = (Power(e1, 2) / ((-e0 + e1) * Power(e1 - e2, 2)) +
              (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * Power(e0 - e2, 2)) -
              (Power(e1, 2) * (e0 * (e1 - 3. * e2) + 2. * e1 * e2) * log_min_e1) /
                  (Power(e0 - e1, 2) * Power(e1 - e2, 3)) +
              (Power(e2, 2) *
               ((e0 - e2) * (e1 - e2) + (3. * e0 * e1 - (e0 + 2. * e1) * e2) * log_min_e2)) /
                  (Power(e0 - e2, 2) * Power(-e1 + e2, 3))) /
             6.;
      r[2] = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
               2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
               (e0 - e1) * e2 *
                   ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                    2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                          (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                        log_min_e2)) /
             (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
      r[3] = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
               2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
               (e0 - e1) * e2 *
                   ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                    2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                          (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                        log_min_e2)) /
             (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_A: {
      assert(not_equal(e0, e1));
      assert(are_equal(e2, e1));
      assert(are_equal(e3, e1));

      r[0] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[1] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
      r[2] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
      r[3] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_B: {
      assert(not_equal(e0, e1));
      assert(are_equal(e0, e2));

      r[0] = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
              6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[1] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[2] = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
              6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[3] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
    } break;

    case math::geometry::FOURFOLD_DEGENERACY: {
      assert(are_equal(e0, e1));
      assert(are_equal(e0, e2));
      assert(are_equal(e0, e3));

      r[0] = 1. / (24. * e0);
      r[1] = 1. / (24. * e0);
      r[2] = 1. / (24. * e0);
      r[3] = 1. / (24. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  std::complex<scalartype> result = 0;

  {
    result += f[0] * r[0];
    result += f[1] * r[1];
    result += f[2] * r[2];
    result += f[3] * r[3];

    assert(result == result);  // make sure there is no NAN
  }

  return result;
}

template <typename scalartype>
std::complex<scalartype> tetrahedron_routines_inverse_matrix_function::integrate_matrix_element_3D(
    math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 4);
  std::complex<scalartype> r[4];
  std::complex<scalartype> f[4];

  int i0 = vec[0].i;
  int i1 = vec[1].i;
  int i2 = vec[2].i;
  int i3 = vec[3].i;

  f[0] = vec[i0].f;
  f[1] = vec[i1].f;
  f[2] = vec[i2].f;
  f[3] = vec[i3].f;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;
  std::complex<scalartype> e2 = vec[2].e;
  std::complex<scalartype> e3 = vec[3].e;

  std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
  std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
  std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
  std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e3));
      assert(not_equal(e3, e0));
      assert(not_equal(e0, e2));
      assert(not_equal(e1, e3));

      r[0] = (-((Power(e0, 2) * (3. * e1 * e2 * e3 + Power(e0, 2) * (e1 + e2 + e3) -
                                 2. * e0 * (e1 * e2 + (e1 + e2) * e3)) *
                 log_min_e0) /
                (Power(e0 - e1, 2) * Power(e0 - e2, 2) * Power(e0 - e3, 2))) +
              (Power(e1, 3) * log_min_e1) / (Power(e0 - e1, 2) * (e1 - e2) * (e1 - e3)) +
              (Power(e2, 3) * log_min_e2) / (Power(e0 - e2, 2) * (-e1 + e2) * (e2 - e3)) +
              ((Power(e0, 2) * (e0 - e3)) / ((e0 - e1) * (e0 - e2)) +
               (Power(e3, 3) * log_min_e3) / ((-e1 + e3) * (-e2 + e3))) /
                  Power(e0 - e3, 2)) /
             6.;
      r[1] = (Power(e1, 2) / ((-e0 + e1) * (e1 - e2) * (e1 - e3)) +
              (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * (e0 - e2) * (e0 - e3)) -
              (Power(e1, 2) * (e0 * (Power(e1, 2) + 3. * e2 * e3 - 2. * e1 * (e2 + e3)) +
                               e1 * (-2. * e2 * e3 + e1 * (e2 + e3))) *
               log_min_e1) /
                  (Power(e0 - e1, 2) * Power(e1 - e2, 2) * Power(e1 - e3, 2)) +
              ((Power(e2, 3) * log_min_e2) / (Power(e1 - e2, 2) * (-e0 + e2)) +
               (Power(e3, 3) * log_min_e3) / ((e0 - e3) * Power(e1 - e3, 2))) /
                  (e2 - e3)) /
             6.;
      r[2] = (Power(e0, 3) * log_min_e0 +
              (-(Power(e1, 3) * Power(e0 - e2, 2) * (e0 - e3) * Power(e2 - e3, 2) * log_min_e1) +
               (e0 - e1) *
                   (Power(e2, 2) * (e0 - e3) * (-e1 + e3) *
                        (e0 * (-2. * e1 * e2 + Power(e2, 2) + 3. * e1 * e3 - 2. * e2 * e3) +
                         e2 * (e1 * e2 - 2. * e1 * e3 + e2 * e3)) *
                        log_min_e2 +
                    (-e0 + e2) * (-e1 + e2) * (Power(e2, 2) * (e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                               (e0 - e2) * (e1 - e2) * Power(e3, 3) * log_min_e3))) /
                  (Power(e1 - e2, 2) * (e1 - e3) * Power(e2 - e3, 2))) /
             (6. * (e0 - e1) * Power(e0 - e2, 2) * (e0 - e3));
      r[3] = (Power(e0, 3) * log_min_e0 +
              (-(Power(e1, 3) * (e0 - e2) * Power(e0 - e3, 2) * Power(e2 - e3, 2) * log_min_e1) +
               (e0 - e1) * (Power(e2, 3) * Power(e0 - e3, 2) * Power(e1 - e3, 2) * log_min_e2 +
                            (e0 - e2) * (-e1 + e2) * Power(e3, 2) *
                                ((e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                 (3. * e0 * e1 * e2 - 2. * (e0 * e1 + (e0 + e1) * e2) * e3 +
                                  (e0 + e1 + e2) * Power(e3, 2)) *
                                     log_min_e3))) /
                  ((e1 - e2) * Power(e1 - e3, 2) * Power(e2 - e3, 2))) /
             (6. * (e0 - e1) * (e0 - e2) * Power(e0 - e3, 2));
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e0));
      assert(are_equal(e2, e3));

      r[0] = ((Power(e0, 2) * (e1 - e2) - e0 * Power(e2, 2) + e1 * Power(e2, 2)) /
                  ((e0 - e1) * Power(e0 - e2, 2) * (e1 - e2)) -
              (Power(e0, 2) * (-3. * e1 * e2 + e0 * (e1 + 2. * e2)) * log_min_e0) /
                  (Power(e0 - e1, 2) * Power(e0 - e2, 3)) +
              ((Power(e1, 3) * log_min_e1) / Power(e0 - e1, 2) +
               (Power(e2, 2) * (-3. * e0 * e1 + 2. * e0 * e2 + e1 * e2) * log_min_e2) /
                   Power(e0 - e2, 3)) /
                  Power(e1 - e2, 2)) /
             6.;
      r[1] = (Power(e1, 2) / ((-e0 + e1) * Power(e1 - e2, 2)) +
              (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * Power(e0 - e2, 2)) -
              (Power(e1, 2) * (e0 * (e1 - 3. * e2) + 2. * e1 * e2) * log_min_e1) /
                  (Power(e0 - e1, 2) * Power(e1 - e2, 3)) +
              (Power(e2, 2) *
               ((e0 - e2) * (e1 - e2) + (3. * e0 * e1 - (e0 + 2. * e1) * e2) * log_min_e2)) /
                  (Power(e0 - e2, 2) * Power(-e1 + e2, 3))) /
             6.;
      r[2] = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
               2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
               (e0 - e1) * e2 *
                   ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                    2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                          (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                        log_min_e2)) /
             (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
      r[3] = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
               2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
               (e0 - e1) * e2 *
                   ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                    2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                          (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                        log_min_e2)) /
             (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_A: {
      assert(not_equal(e0, e1));
      assert(are_equal(e2, e1));
      assert(are_equal(e3, e1));

      r[0] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[1] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
      r[2] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
      r[3] = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
              6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
             (36. * Power(e0 - e1, 4));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_B: {
      assert(not_equal(e0, e1));
      assert(are_equal(e0, e2));
      assert(are_equal(e1, e3));

      r[0] = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
              6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[1] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[2] = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
              6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
             (12. * Power(e0 - e1, 4));
      r[3] = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) + Power(e1, 3) +
              6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
             (12. * Power(e0 - e1, 4));
    } break;

    case math::geometry::FOURFOLD_DEGENERACY: {
      assert(are_equal(e0, e1));
      assert(are_equal(e0, e2));
      assert(are_equal(e0, e3));

      r[0] = 1. / (24. * e0);
      r[1] = 1. / (24. * e0);
      r[2] = 1. / (24. * e0);
      r[3] = 1. / (24. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  std::complex<scalartype> result = 0;

  {
    result += f[0] * r[0];
    result += f[1] * r[1];
    result += f[2] * r[2];
    result += f[3] * r[3];

    assert(result == result);  // make sure there is no NAN
  }

  return result;
}

template <typename scalartype>
void tetrahedron_routines_inverse_matrix_function::integrate_eigenvalues_3D(
    math::geometry::TetrahedronEigenvalueDegeneracy degeneracy,
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 4);

  int i0 = vec[0].i;
  int i1 = vec[1].i;
  int i2 = vec[2].i;
  int i3 = vec[3].i;

  std::complex<scalartype> e0 = vec[0].e;
  std::complex<scalartype> e1 = vec[1].e;
  std::complex<scalartype> e2 = vec[2].e;
  std::complex<scalartype> e3 = vec[3].e;

  std::complex<scalartype> log_min_e0 = vec[0].log_min_e;
  std::complex<scalartype> log_min_e1 = vec[1].log_min_e;
  std::complex<scalartype> log_min_e2 = vec[2].log_min_e;
  std::complex<scalartype> log_min_e3 = vec[3].log_min_e;

  switch (degeneracy) {
    case math::geometry::NO_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e3));
      assert(not_equal(e3, e0));
      assert(not_equal(e0, e2));
      assert(not_equal(e1, e3));

      vec[i0].r = (-((Power(e0, 2) * (3. * e1 * e2 * e3 + Power(e0, 2) * (e1 + e2 + e3) -
                                      2. * e0 * (e1 * e2 + (e1 + e2) * e3)) *
                      log_min_e0) /
                     (Power(e0 - e1, 2) * Power(e0 - e2, 2) * Power(e0 - e3, 2))) +
                   (Power(e1, 3) * log_min_e1) / (Power(e0 - e1, 2) * (e1 - e2) * (e1 - e3)) +
                   (Power(e2, 3) * log_min_e2) / (Power(e0 - e2, 2) * (-e1 + e2) * (e2 - e3)) +
                   ((Power(e0, 2) * (e0 - e3)) / ((e0 - e1) * (e0 - e2)) +
                    (Power(e3, 3) * log_min_e3) / ((-e1 + e3) * (-e2 + e3))) /
                       Power(e0 - e3, 2)) /
                  6.;
      vec[i1].r = (Power(e1, 2) / ((-e0 + e1) * (e1 - e2) * (e1 - e3)) +
                   (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * (e0 - e2) * (e0 - e3)) -
                   (Power(e1, 2) * (e0 * (Power(e1, 2) + 3. * e2 * e3 - 2. * e1 * (e2 + e3)) +
                                    e1 * (-2. * e2 * e3 + e1 * (e2 + e3))) *
                    log_min_e1) /
                       (Power(e0 - e1, 2) * Power(e1 - e2, 2) * Power(e1 - e3, 2)) +
                   ((Power(e2, 3) * log_min_e2) / (Power(e1 - e2, 2) * (-e0 + e2)) +
                    (Power(e3, 3) * log_min_e3) / ((e0 - e3) * Power(e1 - e3, 2))) /
                       (e2 - e3)) /
                  6.;
      vec[i2].r =
          (Power(e0, 3) * log_min_e0 +
           (-(Power(e1, 3) * Power(e0 - e2, 2) * (e0 - e3) * Power(e2 - e3, 2) * log_min_e1) +
            (e0 - e1) *
                (Power(e2, 2) * (e0 - e3) * (-e1 + e3) *
                     (e0 * (-2. * e1 * e2 + Power(e2, 2) + 3. * e1 * e3 - 2. * e2 * e3) +
                      e2 * (e1 * e2 - 2. * e1 * e3 + e2 * e3)) *
                     log_min_e2 +
                 (-e0 + e2) * (-e1 + e2) * (Power(e2, 2) * (e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                            (e0 - e2) * (e1 - e2) * Power(e3, 3) * log_min_e3))) /
               (Power(e1 - e2, 2) * (e1 - e3) * Power(e2 - e3, 2))) /
          (6. * (e0 - e1) * Power(e0 - e2, 2) * (e0 - e3));
      vec[i3].r = (Power(e0, 3) * log_min_e0 +
                   (-(Power(e1, 3) * (e0 - e2) * Power(e0 - e3, 2) * Power(e2 - e3, 2) * log_min_e1) +
                    (e0 - e1) * (Power(e2, 3) * Power(e0 - e3, 2) * Power(e1 - e3, 2) * log_min_e2 +
                                 (e0 - e2) * (-e1 + e2) * Power(e3, 2) *
                                     ((e0 - e3) * (-e1 + e3) * (-e2 + e3) +
                                      (3. * e0 * e1 * e2 - 2. * (e0 * e1 + (e0 + e1) * e2) * e3 +
                                       (e0 + e1 + e2) * Power(e3, 2)) *
                                          log_min_e3))) /
                       ((e1 - e2) * Power(e1 - e3, 2) * Power(e2 - e3, 2))) /
                  (6. * (e0 - e1) * (e0 - e2) * Power(e0 - e3, 2));
    } break;

    case math::geometry::TWOFOLD_DEGENERACY: {
      assert(not_equal(e0, e1));
      assert(not_equal(e1, e2));
      assert(not_equal(e2, e0));
      assert(are_equal(e2, e3));

      vec[i0].r = ((Power(e0, 2) * (e1 - e2) - e0 * Power(e2, 2) + e1 * Power(e2, 2)) /
                       ((e0 - e1) * Power(e0 - e2, 2) * (e1 - e2)) -
                   (Power(e0, 2) * (-3. * e1 * e2 + e0 * (e1 + 2. * e2)) * log_min_e0) /
                       (Power(e0 - e1, 2) * Power(e0 - e2, 3)) +
                   ((Power(e1, 3) * log_min_e1) / Power(e0 - e1, 2) +
                    (Power(e2, 2) * (-3. * e0 * e1 + 2. * e0 * e2 + e1 * e2) * log_min_e2) /
                        Power(e0 - e2, 3)) /
                       Power(e1 - e2, 2)) /
                  6.;
      vec[i1].r = (Power(e1, 2) / ((-e0 + e1) * Power(e1 - e2, 2)) +
                   (Power(e0, 3) * log_min_e0) / (Power(e0 - e1, 2) * Power(e0 - e2, 2)) -
                   (Power(e1, 2) * (e0 * (e1 - 3. * e2) + 2. * e1 * e2) * log_min_e1) /
                       (Power(e0 - e1, 2) * Power(e1 - e2, 3)) +
                   (Power(e2, 2) *
                    ((e0 - e2) * (e1 - e2) + (3. * e0 * e1 - (e0 + 2. * e1) * e2) * log_min_e2)) /
                       (Power(e0 - e2, 2) * Power(-e1 + e2, 3))) /
                  6.;
      vec[i2].r = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
                    2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
                    (e0 - e1) * e2 *
                        ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                         2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                               (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                             log_min_e2)) /
                  (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
      vec[i3].r = -(2. * Power(e0, 3) * Power(e1 - e2, 3) * log_min_e0 -
                    2. * Power(e1, 3) * Power(e0 - e2, 3) * log_min_e1 +
                    (e0 - e1) * e2 *
                        ((e0 - e2) * (e1 - e2) * (5. * e0 * e1 - 3. * (e0 + e1) * e2 + Power(e2, 2)) +
                         2. * (3. * Power(e0, 2) * Power(e1, 2) - 3. * e0 * e1 * (e0 + e1) * e2 +
                               (Power(e0, 2) + e0 * e1 + Power(e1, 2)) * Power(e2, 2)) *
                             log_min_e2)) /
                  (12. * (e0 - e1) * Power(e0 - e2, 3) * Power(-e1 + e2, 3));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_A: {
      assert(not_equal(e0, e1));
      assert(are_equal(e2, e1));
      assert(are_equal(e3, e1));

      vec[i0].r = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) +
                   Power(e1, 3) + 6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
                  (12. * Power(e0 - e1, 4));
      vec[i1].r = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
                   6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
                  (36. * Power(e0 - e1, 4));
      vec[i2].r = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
                   6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
                  (36. * Power(e0 - e1, 4));
      vec[i3].r = (-((e0 - e1) * (11. * Power(e0, 2) - 7. * e0 * e1 + 2. * Power(e1, 2))) +
                   6. * Power(e0, 3) * (log_min_e0 - log_min_e1)) /
                  (36. * Power(e0 - e1, 4));
    } break;

    case math::geometry::THREEFOLD_DEGENERACY_B: {
      assert(not_equal(e0, e1));
      assert(are_equal(e0, e2));
      assert(are_equal(e1, e3));

      vec[i0].r = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
                   6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
                  (12. * Power(e0 - e1, 4));
      vec[i1].r = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) +
                   Power(e1, 3) + 6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
                  (12. * Power(e0 - e1, 4));
      vec[i2].r = ((e0 - e1) * (Power(e0, 2) - 5. * e0 * e1 - 2. * Power(e1, 2)) +
                   6. * e0 * Power(e1, 2) * (log_min_e0 - log_min_e1)) /
                  (12. * Power(e0 - e1, 4));
      vec[i3].r = (2. * Power(e0, 3) + 3. * Power(e0, 2) * e1 - 6. * e0 * Power(e1, 2) +
                   Power(e1, 3) + 6. * Power(e0, 2) * e1 * (-log_min_e0 + log_min_e1)) /
                  (12. * Power(e0 - e1, 4));
    } break;

    case math::geometry::FOURFOLD_DEGENERACY: {
      assert(are_equal(e0, e1));
      assert(are_equal(e0, e2));
      assert(are_equal(e0, e3));

      vec[i0].r = 1. / (24. * e0);
      vec[i1].r = 1. / (24. * e0);
      vec[i2].r = 1. / (24. * e0);
      vec[i3].r = 1. / (24. * e0);
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalartype>
math::geometry::TetrahedronEigenvalueDegeneracy tetrahedron_routines_inverse_matrix_function::find_degeneracy_3D(
    std::vector<matrix_element_struct<scalartype>>& vec) {
  assert(vec.size() == 4);

  if (not_equal(vec[1].e, vec[0].e) and not_equal(vec[2].e, vec[1].e) and
      not_equal(vec[3].e, vec[2].e) and not_equal(vec[0].e, vec[3].e) and
      not_equal(vec[3].e, vec[1].e) and not_equal(vec[2].e, vec[0].e)) {
    return math::geometry::NO_DEGENERACY;
  }

  if (are_equal(vec[1].e, vec[0].e) and are_equal(vec[2].e, vec[0].e) and
      are_equal(vec[3].e, vec[0].e)) {
    return math::geometry::FOURFOLD_DEGENERACY;
  }

  vec[0].i = 0;
  vec[1].i = 1;
  vec[2].i = 2;
  vec[3].i = 3;

  do {
    if (not_equal(vec[1].e, vec[0].e) and not_equal(vec[2].e, vec[1].e) and
        not_equal(vec[0].e, vec[2].e) and are_equal(vec[3].e, vec[2].e)) {
      return math::geometry::TWOFOLD_DEGENERACY;
    }

    if (not_equal(vec[1].e, vec[0].e) and are_equal(vec[2].e, vec[1].e) and
        are_equal(vec[3].e, vec[1].e)) {
      return math::geometry::THREEFOLD_DEGENERACY_A;
    }

    if (not_equal(vec[1].e, vec[0].e) and are_equal(vec[2].e, vec[0].e) and
        are_equal(vec[3].e, vec[1].e)) {
      return math::geometry::THREEFOLD_DEGENERACY_B;
    }
  } while (std::next_permutation(
      vec.begin(), vec.end(),
      tetrahedron_routines_inverse_matrix_function::permutation_comp<scalartype>));

  throw std::logic_error(__FUNCTION__);
  return math::geometry::NO_DEGENERACY;
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_TETRAHEDRON_ROUTINES_INVERSE_MATRIX_FUNCTION_HPP
