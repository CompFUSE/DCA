// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements point_group_symmetry_element.hpp.

#include "dca/phys/domains/quantum/point_group_symmetry_element.hpp"

#include <cmath>
#include <cstring>
#include <numeric>
#include <sstream>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

point_group_symmetry_element::point_group_symmetry_element(int d)
    : DIMENSION(d),

      ORDER(1),
      PHASE(1.),

      O(NULL),
      t(NULL) {
  O = new double[DIMENSION * DIMENSION];
  t = new double[DIMENSION];

  std::memset(t, 0, DIMENSION * sizeof(double));
  std::memset(O, 0, DIMENSION * DIMENSION * sizeof(double));

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      O[i + j * DIMENSION] = i == j ? 1 : 0;
}

point_group_symmetry_element::point_group_symmetry_element(const point_group_symmetry_element& other)
    : DIMENSION(other.DIMENSION),

      ORDER(other.ORDER),
      PHASE(other.PHASE),

      O(NULL),
      t(NULL) {
  P = other.P;

  O = new double[DIMENSION * DIMENSION];
  t = new double[DIMENSION];

  std::memcpy(t, other.t, DIMENSION * sizeof(double));
  std::memcpy(O, other.O, DIMENSION * DIMENSION * sizeof(double));
}

point_group_symmetry_element::~point_group_symmetry_element() {
  delete[] O;
  delete[] t;
}

void point_group_symmetry_element::linear_transform(const double *t0, double* t1) const {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];
}

void point_group_symmetry_element::transform(const double* t0, double* t1) const {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];

  for (int i = 0; i < DIMENSION; ++i)
    t1[i] += t[i];
}

std::string point_group_symmetry_element::describe() const {
  std::ostringstream label;

  if (DIMENSION == 2) {
    // O is column major, so its first column (O[0], O[1]) is (cos, sin) of the rotation angle for a
    // proper rotation (det +1) and (cos, sin) of twice the mirror-line angle for a reflection (det
    // -1). Each operation gets a Schoenflies name (Cn^k for rotation, and sigma for reflection) and
    // a geometric description.
    //
    // Rotations are decoded in twelfths of a full turn, with the parenthetical is still printed in
    // degrees (t * 30) for readability.
    constexpr int kTwelfths = 12;
    const double det = O[0] * O[3] - O[2] * O[1];
    if (det > 0.) {
      const long raw = std::lround(std::atan2(O[1], O[0]) * kTwelfths / (2. * M_PI));
      const int t = static_cast<int>((raw % kTwelfths + kTwelfths) % kTwelfths);
      if (t == 0) {
        label << "E (identity)";
      }
      else {
        const int g = std::gcd(t, kTwelfths);
        label << "C" << kTwelfths / g;
        if (t / g > 1)
          label << "^" << t / g;
        label << " (rotation " << t * (360 / kTwelfths) << " deg)";
      }
    }
    else {
      const long raw = std::lround(std::atan2(O[1], O[0]) / 2. * 180. / M_PI);
      const int deg = static_cast<int>((raw % 180 + 180) % 180);
      label << "sigma (mirror @ " << deg << " deg)";
    }
  }
  else {
    label << DIMENSION << "D symmetry element";
  }

  return label.str();
}

}  // namespace domains
}  // namespace phys
}  // namespace dca
