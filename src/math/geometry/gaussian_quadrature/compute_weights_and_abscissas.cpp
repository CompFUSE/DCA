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
// This file implements the computeWeightsAndAbscissas functions for 1D, 2D and 3D.

#include "dca/math/geometry/gaussian_quadrature/compute_weights_and_abscissas.hpp"

namespace dca {
namespace math {
namespace gaussquad {
// dca::math::gaussquad::

// 1D
template <>
void computeWeightsAndAbscissas(const int rule, math_algorithms::tetrahedron<1>& tet) {
  const int DIM = 1;

  tet.N_q = gm_rule_size(rule, DIM);

  if (tet.q_w != NULL)
    delete[] tet.q_w;
  if (tet.q_vecs != NULL)
    delete[] tet.q_vecs;

  tet.q_w = new double[tet.N_q];
  tet.q_vecs = new double[tet.N_q * DIM];

  double* q_vecs = new double[tet.N_q * DIM];

  gm_rule_set(rule, DIM, tet.N_q, tet.q_w, q_vecs);

  for (int l = 0; l < tet.N_q; l++)
    tet.q_vecs[0 + l * 1] = tet.vec_0[0] + (tet.vec_1[0] - tet.vec_0[0]) * q_vecs[0 + l * DIM];

  delete[] q_vecs;
}

// 2D
template <>
void computeWeightsAndAbscissas(const int rule, math_algorithms::tetrahedron<2>& tet) {
  const int DIM = 2;

  tet.N_q = gm_rule_size(rule, DIM);

  tet.q_w = new double[tet.N_q];
  tet.q_vecs = new double[tet.N_q * DIM];

  double* q_vecs = new double[tet.N_q * DIM];

  gm_rule_set(rule, DIM, tet.N_q, tet.q_w, q_vecs);

  for (int l = 0; l < tet.N_q; l++) {
    tet.q_vecs[0 + l * DIM] = tet.vec_0[0] + (tet.vec_1[0] - tet.vec_0[0]) * q_vecs[0 + l * DIM] +
                              (tet.vec_2[0] - tet.vec_0[0]) * q_vecs[1 + l * DIM];
    tet.q_vecs[1 + l * DIM] = tet.vec_0[1] + (tet.vec_1[1] - tet.vec_0[1]) * q_vecs[0 + l * DIM] +
                              (tet.vec_2[1] - tet.vec_0[1]) * q_vecs[1 + l * DIM];
  }

  delete[] q_vecs;
}

// 3D
template <>
void computeWeightsAndAbscissas(const int rule, math_algorithms::tetrahedron<3>& tet) {
  const int DIM = 3;

  tet.N_q = gm_rule_size(rule, DIM);

  tet.q_w = new double[tet.N_q];
  tet.q_vecs = new double[tet.N_q * DIM];

  double* q_vecs = new double[tet.N_q * DIM];

  gm_rule_set(rule, DIM, tet.N_q, tet.q_w, q_vecs);

  for (int l = 0; l < tet.N_q; l++) {
    tet.q_vecs[0 + l * DIM] = tet.vec_0[0] + (tet.vec_1[0] - tet.vec_0[0]) * q_vecs[0 + l * DIM] +
                              (tet.vec_2[0] - tet.vec_0[0]) * q_vecs[1 + l * DIM] +
                              (tet.vec_3[0] - tet.vec_0[0]) * q_vecs[2 + l * DIM];
    tet.q_vecs[1 + l * DIM] = tet.vec_0[1] + (tet.vec_1[1] - tet.vec_0[1]) * q_vecs[0 + l * DIM] +
                              (tet.vec_2[1] - tet.vec_0[1]) * q_vecs[1 + l * DIM] +
                              (tet.vec_3[1] - tet.vec_0[1]) * q_vecs[2 + l * DIM];
    tet.q_vecs[2 + l * DIM] = tet.vec_0[2] + (tet.vec_1[2] - tet.vec_0[2]) * q_vecs[0 + l * DIM] +
                              (tet.vec_2[2] - tet.vec_0[2]) * q_vecs[1 + l * DIM] +
                              (tet.vec_3[2] - tet.vec_0[2]) * q_vecs[2 + l * DIM];
  }

  delete[] q_vecs;
}

}  // gaussquad
}  // math
}  // dca
