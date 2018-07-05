// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements the tetrahedron class.

#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

//
// Full template specialization for a 1D tetrahedron.
//
tetrahedron<1>::tetrahedron() : N_q(-1), q_w(nullptr), q_vecs(nullptr) {}

tetrahedron<1>::~tetrahedron() {
  if (q_w != nullptr)
    delete[] q_w;

  if (q_vecs != nullptr)
    delete[] q_vecs;
}

std::vector<double> tetrahedron<1>::compute_cm() {
  std::vector<double> cm(1, 0);
  cm[0] = (vec_0[0] + vec_1[0]) / 2.;
  return cm;
}

void tetrahedron<1>::translate(const std::vector<double>& k_point) {
  assert(k_point.size() == 1);

  for (int j = 0; j < 1; j++) {
    vec_0[j] += k_point[j];
    vec_1[j] += k_point[j];
  }

  for (int i = 0; i < N_q; i++)
    for (int j = 0; j < 1; j++)
      q_vecs[j + i * 1] += k_point[j];
}

void tetrahedron<1>::plot(::dca::util::Plot& plot) {
  std::vector<double> x(2);
  std::vector<double> y(2);

  x[0] = vec_0[0];
  y[0] = 0;
  x[1] = vec_1[0];
  y[1] = 0;

  plot.plot(x, y);
}

void tetrahedron<1>::plot_q_vecs(::dca::util::Plot& plot) {
  std::vector<double> x(N_q);
  std::vector<double> y(N_q);

  for (int i = 0; i < N_q; i++) {
    x[i] = q_vecs[i];
    y[i] = 0;
  }

  plot.plot(x, y);
}

//
// Full template specialization for a 2D tetrahedron.
//
tetrahedron<2>::tetrahedron() : N_q(-1), q_w(nullptr), q_vecs(nullptr) {}

tetrahedron<2>::~tetrahedron() {
  if (q_w != nullptr)
    delete[] q_w;

  if (q_vecs != nullptr)
    delete[] q_vecs;
}

std::vector<double> tetrahedron<2>::compute_cm() {
  std::vector<double> cm(2, 0);

  cm[0] = (vec_0[0] + vec_1[0] + vec_2[0]) / 3.;
  cm[1] = (vec_0[1] + vec_1[1] + vec_2[1]) / 3.;

  return cm;
}

void tetrahedron<2>::translate(const std::vector<double>& k_point) {
  assert(k_point.size() == 2);

  for (int j = 0; j < 2; j++) {
    vec_0[j] += k_point[j];
    vec_1[j] += k_point[j];
    vec_2[j] += k_point[j];
  }

  for (int i = 0; i < N_q; i++)
    for (int j = 0; j < 2; j++)
      q_vecs[j + i * 2] += k_point[j];
}

void tetrahedron<2>::plot(::dca::util::Plot& plot) {
  plot.plotLine2D(vec_0, vec_1);
  plot.plotLine2D(vec_1, vec_2);
  plot.plotLine2D(vec_2, vec_0);
}

void tetrahedron<2>::plot_q_vecs(::dca::util::Plot& plot) {
  if (q_vecs != nullptr) {
    std::vector<double> x(N_q);
    std::vector<double> y(N_q);

    for (int i = 0; i < N_q; i++) {
      x[i] = q_vecs[0 + i * 2];
      y[i] = q_vecs[1 + i * 2];
    }

    plot.plot(x, y);
  }
}

bool tetrahedron<2>::pair_same(std::pair<std::complex<double>, std::complex<double>> const& x,
                               std::pair<std::complex<double>, std::complex<double>> const& y) {
  double abs_x = abs(x.first);
  double abs_y = abs(y.first);

  if (abs_x < 1. && abs_y < 1.)
    return abs(x.first - y.first) < 1.e-6;
  else {
    double MAX = abs_x > abs_y ? abs_x : abs_y;
    return abs(x.first - y.first) < ((1.e-6) * MAX);
  }
}

bool tetrahedron<2>::pair_same_index(std::pair<std::complex<double>, int> const& x,
                                     std::pair<std::complex<double>, int> const& y) {
  double abs_x = abs(x.first);
  double abs_y = abs(y.first);

  if (abs_x < 1. && abs_y < 1.)
    return abs(x.first - y.first) < 1.e-6;
  else {
    double MAX = abs_x > abs_y ? abs_x : abs_y;
    return abs(x.first - y.first) < ((1.e-6) * MAX);
  }
}

//
// Full template specialization for a 3D tetrahedron.
//
tetrahedron<3>::tetrahedron() : N_q(-1), q_w(nullptr), q_vecs(nullptr) {}

tetrahedron<3>::~tetrahedron() {
  if (q_w != nullptr)
    delete[] q_w;

  if (q_vecs != nullptr)
    delete[] q_vecs;
}

std::vector<double> tetrahedron<3>::compute_cm() {
  std::vector<double> cm(3, 0);

  cm[0] = (vec_0[0] + vec_1[0] + vec_2[0] + vec_3[0]) / 4.;
  cm[1] = (vec_0[1] + vec_1[1] + vec_2[1] + vec_3[1]) / 4.;
  cm[2] = (vec_0[2] + vec_1[2] + vec_2[2] + vec_3[2]) / 4.;

  return cm;
}

void tetrahedron<3>::translate(const std::vector<double>& k_point) {
  assert(k_point.size() == 3);

  for (int j = 0; j < 3; j++) {
    vec_0[j] += k_point[j];
    vec_1[j] += k_point[j];
    vec_2[j] += k_point[j];
    vec_3[j] += k_point[j];
  }

  for (int i = 0; i < N_q; i++)
    for (int j = 0; j < 3; j++)
      q_vecs[j + i * 3] += k_point[j];
}

void tetrahedron<3>::plot(::dca::util::Plot& plot) {
  plot.plotLine3D(vec_0, vec_1);
  plot.plotLine3D(vec_1, vec_2);
  plot.plotLine3D(vec_2, vec_3);
  plot.plotLine3D(vec_3, vec_0);
  plot.plotLine3D(vec_0, vec_2);
  plot.plotLine3D(vec_1, vec_3);
}

bool tetrahedron<3>::pair_same(std::pair<std::complex<double>, std::complex<double>> const& x,
                               std::pair<std::complex<double>, std::complex<double>> const& y) {
  double abs_x = norm(x.first);
  double abs_y = norm(y.first);

  if (abs_x < 1. && abs_y < 1.)
    return norm(x.first - y.first) < 1.e-3;
  else {
    double MAX = abs_x > abs_y ? abs_x : abs_y;
    return norm(x.first - y.first) < ((1.e-3) * MAX);
  }
}

}  // geometry
}  // math
}  // dca
