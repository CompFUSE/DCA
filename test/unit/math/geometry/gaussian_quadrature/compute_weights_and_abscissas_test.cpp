// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the computeWeightsAndAbscissas functions for the Gaussian quadrature algorithm.

#include "dca/math/geometry/gaussian_quadrature/compute_weights_and_abscissas.hpp"
#include "gtest/gtest.h"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"

TEST(ComputeWeightsAndAbscissasTest, 1DUnitSimplex) {
  // 1D unit simplex = line segment of length 1
  dca::math::geometry::tetrahedron<1> tet;
  tet.vec_0 = {0.};
  tet.vec_1 = {1.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(3, tet.N_q);

  EXPECT_FLOAT_EQ(0.66666669, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.66666669, tet.q_w[1]);
  EXPECT_FLOAT_EQ(-0.33333334, tet.q_w[2]);

  EXPECT_FLOAT_EQ(0.25, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(0.75, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(0.5, tet.q_vecs[2]);
}

TEST(ComputeWeightsAndAbscissasTest, 1DTransformedSimplex) {
  // Transformed 1D simplex: unit simplex streched by a factor of 2 and shifted by 1
  dca::math::geometry::tetrahedron<1> tet;
  tet.vec_0 = {1.};
  tet.vec_1 = {3.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(3, tet.N_q);

  EXPECT_FLOAT_EQ(0.66666669, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.66666669, tet.q_w[1]);
  EXPECT_FLOAT_EQ(-0.33333334, tet.q_w[2]);

  EXPECT_FLOAT_EQ(1.5, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(2.5, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[2]);
}

TEST(ComputeWeightsAndAbscissasTest, 2DUnitSimplex) {
  // 2D unit simplex = isosceles triangle of side length 1
  dca::math::geometry::tetrahedron<2> tet;
  tet.vec_0 = {0., 0.};
  tet.vec_1 = {1., 0.};
  tet.vec_2 = {0., 1.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(4, tet.N_q);

  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[1]);
  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[2]);
  EXPECT_FLOAT_EQ(-0.5625, tet.q_w[3]);

  EXPECT_FLOAT_EQ(0.2, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(0.2, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(0.6, tet.q_vecs[2]);
  EXPECT_FLOAT_EQ(0.2, tet.q_vecs[3]);
  EXPECT_FLOAT_EQ(0.2, tet.q_vecs[4]);
  EXPECT_FLOAT_EQ(0.6, tet.q_vecs[5]);
  EXPECT_FLOAT_EQ(0.33333334, tet.q_vecs[6]);
  EXPECT_FLOAT_EQ(0.33333334, tet.q_vecs[7]);
}

TEST(ComputeWeightsAndAbscissasTest, 2DTransformedSimplex) {
  dca::math::geometry::tetrahedron<2> tet;
  tet.vec_0 = {2., 1.};
  tet.vec_1 = {3., 4.};
  tet.vec_2 = {0., 3.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(4, tet.N_q);

  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[1]);
  EXPECT_FLOAT_EQ(0.52083331, tet.q_w[2]);
  EXPECT_FLOAT_EQ(-0.5625, tet.q_w[3]);

  EXPECT_FLOAT_EQ(1.8, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(2.2, tet.q_vecs[2]);
  EXPECT_FLOAT_EQ(3.2, tet.q_vecs[3]);
  EXPECT_FLOAT_EQ(1.0, tet.q_vecs[4]);
  EXPECT_FLOAT_EQ(2.8, tet.q_vecs[5]);
  EXPECT_FLOAT_EQ(1.6666666, tet.q_vecs[6]);
  EXPECT_FLOAT_EQ(2.6666667, tet.q_vecs[7]);
}

TEST(ComputeWeightsAndAbscissasTest, 3DUnitSimplex) {
  dca::math::geometry::tetrahedron<3> tet;
  tet.vec_0 = {0., 0., 0.};
  tet.vec_1 = {1., 0., 0.};
  tet.vec_2 = {0., 1., 0.};
  tet.vec_3 = {0., 0., 1.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(5, tet.N_q);

  EXPECT_FLOAT_EQ(0.45, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[1]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[2]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[3]);
  EXPECT_FLOAT_EQ(-0.8, tet.q_w[4]);

  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[2]);
  EXPECT_FLOAT_EQ(0.5, tet.q_vecs[3]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[4]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[5]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[6]);
  EXPECT_FLOAT_EQ(0.5, tet.q_vecs[7]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[8]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[9]);
  EXPECT_FLOAT_EQ(0.16666667, tet.q_vecs[10]);
  EXPECT_FLOAT_EQ(0.5, tet.q_vecs[11]);
  EXPECT_FLOAT_EQ(0.25, tet.q_vecs[12]);
  EXPECT_FLOAT_EQ(0.25, tet.q_vecs[13]);
  EXPECT_FLOAT_EQ(0.25, tet.q_vecs[14]);
}

TEST(ComputeWeightsAndAbscissasTest, 3DTransformedSimplex) {
  dca::math::geometry::tetrahedron<3> tet;
  tet.vec_0 = {2., 1., 1.};
  tet.vec_1 = {3., 4., 0.};
  tet.vec_2 = {0., 4., 2.};
  tet.vec_3 = {0., 1., 5.};

  const int rule = 1;

  dca::math::geometry::computeWeightsAndAbscissas(rule, tet);

  EXPECT_EQ(5, tet.N_q);

  EXPECT_FLOAT_EQ(0.45, tet.q_w[0]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[1]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[2]);
  EXPECT_FLOAT_EQ(0.45, tet.q_w[3]);
  EXPECT_FLOAT_EQ(-0.8, tet.q_w[4]);

  EXPECT_FLOAT_EQ(1.5, tet.q_vecs[0]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[1]);
  EXPECT_FLOAT_EQ(1.6666666, tet.q_vecs[2]);
  EXPECT_FLOAT_EQ(1.8333334, tet.q_vecs[3]);
  EXPECT_FLOAT_EQ(3.0, tet.q_vecs[4]);
  EXPECT_FLOAT_EQ(1.3333334, tet.q_vecs[5]);
  EXPECT_FLOAT_EQ(0.83333331, tet.q_vecs[6]);
  EXPECT_FLOAT_EQ(3.0, tet.q_vecs[7]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[8]);
  EXPECT_FLOAT_EQ(0.83333331, tet.q_vecs[9]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[10]);
  EXPECT_FLOAT_EQ(3.0, tet.q_vecs[11]);
  EXPECT_FLOAT_EQ(1.25, tet.q_vecs[12]);
  EXPECT_FLOAT_EQ(2.5, tet.q_vecs[13]);
  EXPECT_FLOAT_EQ(2.0, tet.q_vecs[14]);
}
