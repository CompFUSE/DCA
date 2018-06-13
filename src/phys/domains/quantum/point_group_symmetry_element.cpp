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

#include <cstring>

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

void point_group_symmetry_element::linear_transform(double* t0, double* t1) {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];
}

void point_group_symmetry_element::transform(double* t0, double* t1) {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];

  for (int i = 0; i < DIMENSION; ++i)
    t1[i] += t[i];
}

}  // domains
}  // phys
}  // dca
