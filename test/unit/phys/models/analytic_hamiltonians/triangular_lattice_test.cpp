// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests triangular_lattice.hpp.

#include "dca/phys/models/analytic_hamiltonians/triangular_lattice.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_hexagonal.hpp"
#include "dca/phys/parameters/model_parameters.hpp"

TEST(TriangularLatticeTest, Initialize_H_0) {
  using PointGroup = dca::phys::domains::C6;
  using Lattice = dca::phys::models::triangular_lattice<PointGroup>;

  using BandDmn = dca::func::dmn<1, int>;
  using SpinDmn = dca::func::dmn<2, int>;
  using BandSpinDmn = dca::func::dmn_variadic<dca::func::dmn_0<BandDmn>, dca::func::dmn_0<SpinDmn>>;

  using KDmn = dca::func::dmn<4, std::vector<double>>;
  KDmn::set_elements({{0., 0.},
                      {M_PI, -std::sin(M_PI / 3.) * 2. / 3. * M_PI},
                      {M_PI, std::sin(M_PI / 3.) * 2. / 3. * M_PI},
                      {2. * M_PI, 0.}});

  dca::func::function<std::complex<double>,
                      dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, dca::func::dmn_0<KDmn>>>
      H_0;

  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<Lattice>> params;
  params.set_t(1.);

  Lattice::initialize_H_0(params, H_0);

  // All imaginary parts should be zero.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int s1 = 0; s1 < SpinDmn::dmn_size(); ++s1)
      for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
        for (int s2 = 0; s2 < SpinDmn::dmn_size(); ++s2)
          for (int k = 0; k < KDmn::dmn_size(); ++k)
            EXPECT_DOUBLE_EQ(0., H_0(b1, s1, b2, s2, k).imag());

  // All matrix elements with different spin indices should be zero.
  for (int b1 = 0; b1 < BandDmn::dmn_size(); ++b1)
    for (int b2 = 0; b2 < BandDmn::dmn_size(); ++b2)
      for (int k = 0; k < KDmn::dmn_size(); ++k) {
        EXPECT_DOUBLE_EQ(0., H_0(b1, 0, b2, 1, k).real());
        EXPECT_DOUBLE_EQ(0., H_0(b1, 1, b2, 0, k).real());
      }

  // Check nonvanishing Hamiltonian matrix elements.
  EXPECT_DOUBLE_EQ(-6., H_0(0, 0, 0, 0, 0).real());
  EXPECT_DOUBLE_EQ(-6., H_0(0, 1, 0, 1, 0).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 1).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 2).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 0, 0, 0, 3).real());
  EXPECT_DOUBLE_EQ(2., H_0(0, 1, 0, 1, 3).real());
}
