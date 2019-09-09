// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements tetrahedron_routines_harmonic_function.hpp.

#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_routines_harmonic_function.hpp"

#include <cassert>
#include <cmath>

#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

/************************************
 ***
 ***   1D harmonic-integration
 ***
 ************************************/

std::complex<double> tetrahedron_routines_harmonic_function::execute(
    const std::vector<double>& r_vec, const math::geometry::tetrahedron<1>& tetrahedron) {
  assert(r_vec.size() == 1);

  double r = r_vec[0];

  std::vector<double> K0 = tetrahedron.vec_0;
  std::vector<double> K1 = tetrahedron.vec_0;

  double a0 = K0[0];
  double a1 = K1[0];

  std::complex<double> result(0., 0.);

  if (std::abs(r) < 1.e-6) {
    result.real(a1 - a0);
    result.imag(0);
  }
  else {
    result.real(-std::sin(a0 * r) / r + std::sin(a1 * r) / r);
    result.imag(std::cos(a0 * r) / r - std::cos(a1 * r) / r);
  }

  return result;
}

/************************************
 ***
 ***   2D harmonic-integration
 ***
 ************************************/

/*!
 *   We want to integrate a harmonic function e^{i k r} over an arbitrary triangle, defined by k_0,
 * k_1, k_2.
 *   This comes down to:
 *
 *   \int_0^1 dt_1 \int_0^{1-t_1} dt_2 e^{i r k} * det(d_1, d_2)
 *
 *   with
 *          k   = k_0 + d_1*t_1 + d_2*t_2
 *
 *          d_1 = k_1-k_0,
 *          d_2 = k_2-k_0.
 *          d_3 = d_2-d_1 = k_2-k_1
 *
 *   for convinience, we define
 *
 *          r_dot_di = r \dot d_i
 *
 *   there are 4 cases:
 *
 *    case 0:  |r_dot_d1|>0 and |r_dot_d2|>0  and |r_dot_d3|>0 (general)
 *    case 1:  |r_dot_d1|=0 and |r_dot_d2|>0  and |r_dot_d3|>0 (special 1)
 *    case 2:  |r_dot_d1|>0 and |r_dot_d2|=0  and |r_dot_d3|>0 (special 2)
 *    case 3:  |r_dot_d1|>0 and |r_dot_d2|>0  and |r_dot_d3|=0 (special 3)
 *
 *
 */
std::complex<double> tetrahedron_routines_harmonic_function::execute(
    const std::vector<double>& r_vec, const math::geometry::tetrahedron<2>& tetrahedron) {
  assert(r_vec.size() == 2);

  const static std::complex<double> I(0, 1);
  const static double EPSILON = 1.e-6;

  std::complex<double> result(0., 0.);

  if (std::abs(math::util::innerProduct(r_vec, r_vec)) < EPSILON) {
    result.real(tetrahedron.volume);
    result.imag(0);
  }
  else {
    std::vector<double> K0 = tetrahedron.vec_0;

    std::vector<double> D1 = K0;
    std::vector<double> D2 = K0;
    std::vector<double> D2_min_D1 = K0;

    for (int d = 0; d < 2; d++) {
      D1[d] = tetrahedron.vec_1[d] - tetrahedron.vec_0[d];
      D2[d] = tetrahedron.vec_2[d] - tetrahedron.vec_0[d];

      D2_min_D1[d] = D2[d] - D1[d];
    }

    double dot_R_K0 = math::util::innerProduct(r_vec, K0);
    double dot_R_D1 = math::util::innerProduct(r_vec, D1);
    double dot_R_D2 = math::util::innerProduct(r_vec, D2);

    double dot_R_D2_min_D1 = math::util::innerProduct(r_vec, D2_min_D1);

    double det = math::util::area(D1, D2);
    std::complex<double> phase = std::cos(dot_R_K0) + I * std::sin(dot_R_K0);

    if (std::abs(dot_R_D2_min_D1) > EPSILON) {
      if (std::abs(dot_R_D1) > EPSILON and std::abs(dot_R_D2) > EPSILON)
        result = case_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1) * phase * det;

      if (std::abs(dot_R_D1) < EPSILON)
        result = case_d1_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1) * phase * det;

      if (std::abs(dot_R_D2) < EPSILON)
        result = case_d2_2D(dot_R_D1, dot_R_D2, dot_R_D2_min_D1) * phase * det;
    }
    else
    //        if(abs(dot_R_D2_min_D1)<EPSILON)
    {
      math::geometry::tetrahedron<2> tetrahedron_new;

      permute(tetrahedron_new, tetrahedron);

      result = execute(r_vec, tetrahedron_new);
      // result = case_3_2D(dot_R_D1, dot_R_D2, dot_R_D3)*phase*det;
    }
  }

  return result;
}

void tetrahedron_routines_harmonic_function::permute(math::geometry::tetrahedron<2>& tetrahedron_new,
                                                     const math::geometry::tetrahedron<2> &tetrahedron_old) {
  tetrahedron_new.vec_1 = tetrahedron_old.vec_0;
  tetrahedron_new.vec_2 = tetrahedron_old.vec_1;
  tetrahedron_new.vec_0 = tetrahedron_old.vec_2;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_2D(double dotRD1, double dotRD2,
                                                                     double dotRD2minD1) {
  assert(std::abs(dotRD1) > 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  assert(std::abs(dotRD2minD1) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(-(1 / (dotRD1 * dotRD2minD1)) + 1 / (dotRD2 * dotRD2minD1) +
              std::cos(dotRD1) / (dotRD1 * dotRD2minD1) - std::cos(dotRD2) / (dotRD2 * dotRD2minD1));
  result.imag(std::sin(dotRD1) / (dotRD1 * dotRD2minD1) - std::sin(dotRD2) / (dotRD2 * dotRD2minD1));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d1_2D(double /*dotRD1*/,
                                                                        double dotRD2,
                                                                        double dotRD2minD1) {
  // assert(std::abs(dotRD1) < 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  assert(std::abs(dotRD2minD1) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(1 / (dotRD2 * dotRD2minD1) - std::cos(dotRD2) / (dotRD2 * dotRD2minD1));
  result.imag(1 / dotRD2minD1 - std::sin(dotRD2) / (dotRD2 * dotRD2minD1));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d2_2D(double dotRD1,
                                                                        double /*dotRD2*/,
                                                                        double dotRD2minD1) {
  assert(std::abs(dotRD1) > 1.e-6);
  // assert(std::abs(dotRD2) < 1.e-6);
  assert(std::abs(dotRD2minD1) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(-(1 / (dotRD1 * dotRD2minD1)) + std::cos(dotRD1) / (dotRD1 * dotRD2minD1));
  result.imag(-(1 / dotRD2minD1) + std::sin(dotRD1) / (dotRD1 * dotRD2minD1));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

/*
  std::complex<double> tetrahedron_routines_harmonic_function::case_3_2D(double dotRD1,
  double dotRD2,
  double dotRD3)
  {
  assert(std::abs(dotRD1) > 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  assert(std::abs(dotRD3) < 1.e-6);

  std::complex<double> result(0., 0.);

  real(result) = -std::pow(dotRD2,-2) + std::cos(dotRD2)/std::pow(dotRD2,2) +
  std::sin(dotRD2)/dotRD2;
  imag(result) = -(std::cos(dotRD2)/dotRD2) + std::sin(dotRD2)/std::pow(dotRD2,2);

  return result;
  }
*/

/************************************
 ***
 ***   3D harmonic-integration
 ***
 ************************************/

void tetrahedron_routines_harmonic_function::permute(math::geometry::tetrahedron<3>& tetrahedron_new,
                                                     const math::geometry::tetrahedron<3> &tetrahedron_old) {
  tetrahedron_new.vec_1 = tetrahedron_old.vec_0;
  tetrahedron_new.vec_2 = tetrahedron_old.vec_1;
  tetrahedron_new.vec_3 = tetrahedron_old.vec_2;
  tetrahedron_new.vec_0 = tetrahedron_old.vec_3;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_3D(double dotRD1, double dotRD2,
                                                                     double dotRD3,
                                                                     double dotRD2minD1,
                                                                     double dotRD3minD2,
                                                                     double dotRD1minD3) {
  assert(std::abs(dotRD1) > 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  assert(std::abs(dotRD3) > 1.e-6);

  assert(std::abs(dotRD2minD1) > 1.e-6);
  assert(std::abs(dotRD3minD2) > 1.e-6);
  assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(-((dotRD2 * std::sin(dotRD1)) / (dotRD1 * dotRD1minD3 * dotRD2minD1 * dotRD3minD2)) +
              (dotRD3 * std::sin(dotRD1)) / (dotRD1 * dotRD1minD3 * dotRD2minD1 * dotRD3minD2) +
              std::sin(dotRD2) / (dotRD2 * dotRD2minD1 * dotRD3minD2) +
              std::sin(dotRD3) / (dotRD1minD3 * dotRD3 * dotRD3minD2));

  result.imag(-(1 / (dotRD1 * dotRD2 * dotRD3minD2)) + 1 / (dotRD1 * dotRD3 * dotRD3minD2) +
              (dotRD2 * std::cos(dotRD1)) / (dotRD1 * dotRD1minD3 * dotRD2minD1 * dotRD3minD2) -
              (dotRD3 * std::cos(dotRD1)) / (dotRD1 * dotRD1minD3 * dotRD2minD1 * dotRD3minD2) -
              std::cos(dotRD2) / (dotRD2 * dotRD2minD1 * dotRD3minD2) -
              std::cos(dotRD3) / (dotRD1minD3 * dotRD3 * dotRD3minD2));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d1_3D(double /*dotRD1*/,
                                                                        double dotRD2, double dotRD3,
                                                                        double /*dotRD2minD1*/,
                                                                        double dotRD3minD2,
                                                                        double /*dotRD1minD3*/) {
  // assert(std::abs(dotRD1) < 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  assert(std::abs(dotRD3) > 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  if (std::abs(dotRD3minD2) > 1.e-6) {
    result.real(-(1 / (dotRD2 * dotRD3minD2)) + 1 / (dotRD3 * dotRD3minD2) +
                std::sin(dotRD2) / (std::pow(dotRD2, 2) * dotRD3minD2) -
                std::sin(dotRD3) / (std::pow(dotRD3, 2) * dotRD3minD2));

    result.imag(1 / (std::pow(dotRD2, 2) * dotRD3minD2) - 1 / (std::pow(dotRD3, 2) * dotRD3minD2) -
                std::cos(dotRD2) / (std::pow(dotRD2, 2) * dotRD3minD2) +
                std::cos(dotRD3) / (std::pow(dotRD3, 2) * dotRD3minD2));
  }
  else {
    // cout << __FUNCTION__ << " needs implementation\n";

    result.real(-std::pow(dotRD2, -2) - std::cos(dotRD2) / std::pow(dotRD2, 2) +
                (2 * std::sin(dotRD2)) / std::pow(dotRD2, 3));

    result.imag(2 / std::pow(dotRD2, 3) - (2 * std::cos(dotRD2)) / std::pow(dotRD2, 3) -
                std::sin(dotRD2) / std::pow(dotRD2, 2));
  }

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d2_3D(
    double dotRD1, double /*dotRD2*/, double dotRD3, double /*dotRD2minD1*/, double /*dotRD3minD2*/,
    double dotRD1minD3) {
  assert(std::abs(dotRD1) > 1.e-6);
  // assert(std::abs(dotRD2) < 1.e-6);
  assert(std::abs(dotRD3) > 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  if (std::abs(dotRD1minD3) > 1.e-6) {
    result.real(1 / (dotRD1 * dotRD1minD3) - 1 / (dotRD1minD3 * dotRD3) -
                std::sin(dotRD1) / (std::pow(dotRD1, 2) * dotRD1minD3) +
                std::sin(dotRD3) / (dotRD1minD3 * std::pow(dotRD3, 2)));

    result.imag(-(1 / (std::pow(dotRD1, 2) * dotRD1minD3)) + 1 / (dotRD1minD3 * std::pow(dotRD3, 2)) +
                std::cos(dotRD1) / (std::pow(dotRD1, 2) * dotRD1minD3) -
                std::cos(dotRD3) / (dotRD1minD3 * std::pow(dotRD3, 2)));
  }
  else {
    // cout << __FUNCTION__ << " needs implementation\n";

    result.real(-std::pow(dotRD3, -2) - std::cos(dotRD3) / std::pow(dotRD3, 2) +
                (2 * std::sin(dotRD3)) / std::pow(dotRD3, 3));

    result.imag(2 / std::pow(dotRD3, 3) - (2 * std::cos(dotRD3)) / std::pow(dotRD3, 3) -
                std::sin(dotRD3) / std::pow(dotRD3, 2));
  }

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d3_3D(double dotRD1, double dotRD2,
                                                                        double /*dotRD3*/,
                                                                        double dotRD2minD1,
                                                                        double /*dotRD3minD2*/,
                                                                        double /*dotRD1minD3*/) {
  assert(std::abs(dotRD1) > 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  // assert(std::abs(dotRD3) < 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  if (std::abs(dotRD2minD1) > 1.e-6) {
    result.real(-(1 / (dotRD1 * dotRD2minD1)) + 1 / (dotRD2 * dotRD2minD1) +
                std::sin(dotRD1) / (std::pow(dotRD1, 2) * dotRD2minD1) -
                std::sin(dotRD2) / (std::pow(dotRD2, 2) * dotRD2minD1));

    result.imag(1 / (std::pow(dotRD1, 2) * dotRD2minD1) - 1 / (std::pow(dotRD2, 2) * dotRD2minD1) -
                std::cos(dotRD1) / (std::pow(dotRD1, 2) * dotRD2minD1) +
                std::cos(dotRD2) / (std::pow(dotRD2, 2) * dotRD2minD1));
  }
  else {
    // cout << __FUNCTION__ << " needs implementation\n";

    result.real(-std::pow(dotRD1, -2) - std::cos(dotRD1) / std::pow(dotRD1, 2) +
                (2 * std::sin(dotRD1)) / std::pow(dotRD1, 3));

    result.imag(2 / std::pow(dotRD1, 3) - (2 * std::cos(dotRD1)) / std::pow(dotRD1, 3) -
                std::sin(dotRD1) / std::pow(dotRD1, 2));
  }

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d1_d2_3D(
    double /*dotRD1*/, double /*dotRD2*/, double dotRD3, double /*dotRD2minD1*/,
    double /*dotRD3minD2*/, double /*dotRD1minD3*/) {
  // assert(std::abs(dotRD1) < 1.e-6);
  // assert(std::abs(dotRD2) < 1.e-6);
  assert(std::abs(dotRD3) > 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(std::pow(dotRD3, -2) - std::sin(dotRD3) / std::pow(dotRD3, 3));

  result.imag(-std::pow(dotRD3, -3) + 1 / (2. * dotRD3) + std::cos(dotRD3) / std::pow(dotRD3, 3));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d2_d3_3D(
    double dotRD1, double /*dotRD2*/, double /*dotRD3*/, double /*dotRD2minD1*/,
    double /*dotRD3minD2*/, double /*dotRD1minD3*/) {
  assert(std::abs(dotRD1) > 1.e-6);
  // assert(std::abs(dotRD2) < 1.e-6);
  // assert(std::abs(dotRD3) < 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(std::pow(dotRD1, -2) - std::sin(dotRD1) / std::pow(dotRD1, 3));

  result.imag(-std::pow(dotRD1, -3) + 1 / (2. * dotRD1) + std::cos(dotRD1) / std::pow(dotRD1, 3));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::case_d3_d1_3D(
    double /*dotRD1*/, double dotRD2, double /*dotRD3*/, double /*dotRD2minD1*/,
    double /*dotRD3minD2*/, double /*dotRD1minD3*/) {
  // assert(std::abs(dotRD1) < 1.e-6);
  assert(std::abs(dotRD2) > 1.e-6);
  // assert(std::abs(dotRD3) < 1.e-6);

  //     assert(std::abs(dotRD2minD1) > 1.e-6);
  //     assert(std::abs(dotRD3minD2) > 1.e-6);
  //     assert(std::abs(dotRD1minD3) > 1.e-6);

  std::complex<double> result(0., 0.);

  result.real(std::pow(dotRD2, -2) - std::sin(dotRD2) / std::pow(dotRD2, 3));

  result.imag(-std::pow(dotRD2, -3) + 1 / (2. * dotRD2) + std::cos(dotRD2) / std::pow(dotRD2, 3));

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

std::complex<double> tetrahedron_routines_harmonic_function::execute(
    const std::vector<double>& r_vec, const math::geometry::tetrahedron<3>& tetrahedron) {
  assert(r_vec.size() == 3);

  const static std::complex<double> I(0, 1);
  const static double EPSILON = 1.e-6;

  std::complex<double> result(0., 0.);

  if (std::abs(math::util::innerProduct(r_vec, r_vec)) < EPSILON) {
    result.real(tetrahedron.volume);
    result.imag(0);
  }
  else {
    std::vector<double> K_0 = tetrahedron.vec_0;

    std::vector<double> D1 = K_0;
    std::vector<double> D2 = K_0;
    std::vector<double> D3 = K_0;

    std::vector<double> D2_min_D1 = K_0;
    std::vector<double> D3_min_D2 = K_0;
    std::vector<double> D1_min_D3 = K_0;

    for (int d = 0; d < 3; d++) {
      D1[d] = tetrahedron.vec_1[d] - tetrahedron.vec_0[d];
      D2[d] = tetrahedron.vec_2[d] - tetrahedron.vec_0[d];
      D3[d] = tetrahedron.vec_3[d] - tetrahedron.vec_0[d];

      D2_min_D1[d] = D2[d] - D1[d];
      D3_min_D2[d] = D3[d] - D2[d];
      D1_min_D3[d] = D1[d] - D3[d];
    }

    double dot_R_K0 = math::util::innerProduct(r_vec, K_0);

    double dot_R_D1 = math::util::innerProduct(r_vec, D1);
    double dot_R_D2 = math::util::innerProduct(r_vec, D2);
    double dot_R_D3 = math::util::innerProduct(r_vec, D3);

    double dot_R_D2_min_D1 = math::util::innerProduct(r_vec, D2_min_D1);
    double dot_R_D3_min_D2 = math::util::innerProduct(r_vec, D3_min_D2);
    double dot_R_D1_min_D3 = math::util::innerProduct(r_vec, D1_min_D3);

    double det = math::util::volume(D1, D2, D3);
    std::complex<double> phase = std::cos(dot_R_K0) + I * std::sin(dot_R_K0);

    if (std::abs(dot_R_D1) > EPSILON and std::abs(dot_R_D2) > EPSILON and
        std::abs(dot_R_D3) > EPSILON and std::abs(dot_R_D2_min_D1) > EPSILON and
        std::abs(dot_R_D3_min_D2) > EPSILON and std::abs(dot_R_D1_min_D3) > EPSILON)  // general case
    {
      result =
          case_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2, dot_R_D1_min_D3) *
          det * phase;
    }
    else {
      if (std::abs(dot_R_D1) < EPSILON or std::abs(dot_R_D2) < EPSILON or
          std::abs(dot_R_D3) < EPSILON)  // special cases where one or two dot-products are zero
      {
        if (std::abs(dot_R_D1) < EPSILON and std::abs(dot_R_D2) > EPSILON and
            std::abs(dot_R_D3) > EPSILON)
          result = case_d1_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                              dot_R_D1_min_D3) *
                   det * phase;

        if (std::abs(dot_R_D1) > EPSILON and std::abs(dot_R_D2) < EPSILON and
            std::abs(dot_R_D3) > EPSILON)
          result = case_d2_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                              dot_R_D1_min_D3) *
                   det * phase;

        if (std::abs(dot_R_D1) > EPSILON and std::abs(dot_R_D2) > EPSILON and
            std::abs(dot_R_D3) < EPSILON)
          result = case_d3_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                              dot_R_D1_min_D3) *
                   det * phase;

        if (std::abs(dot_R_D1) < EPSILON and std::abs(dot_R_D2) < EPSILON and
            std::abs(dot_R_D3) > EPSILON)
          result = case_d1_d2_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                                 dot_R_D1_min_D3) *
                   det * phase;

        if (std::abs(dot_R_D1) > EPSILON and std::abs(dot_R_D2) < EPSILON and
            std::abs(dot_R_D3) < EPSILON)
          result = case_d2_d3_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                                 dot_R_D1_min_D3) *
                   det * phase;

        if (std::abs(dot_R_D1) < EPSILON and std::abs(dot_R_D2) > EPSILON and
            std::abs(dot_R_D3) < EPSILON)
          result = case_d3_d1_3D(dot_R_D1, dot_R_D2, dot_R_D3, dot_R_D2_min_D1, dot_R_D3_min_D2,
                                 dot_R_D1_min_D3) *
                   det * phase;
      }
      else {
        //                 cout << "\n\t start permuting\t";
        //                 math::util::print(r_vec);
        //                 cout << "\n";
        //                 math::util::print(tetrahedron.vec_0);cout << "\n";
        //                 math::util::print(tetrahedron.vec_1);cout << "\n";
        //                 math::util::print(tetrahedron.vec_2);cout << "\n";
        //                 math::util::print(tetrahedron.vec_3);cout << "\n";

        math::geometry::tetrahedron<3> tetrahedron_new;

        permute(tetrahedron_new, tetrahedron);

        result = execute(r_vec, tetrahedron_new);
      }
    }
  }

  assert(real(result) == real(result));
  assert(imag(result) == imag(result));

  return result;
}

}  // namespace clustermapping
}  // namespace phys
}  // namespace dca
