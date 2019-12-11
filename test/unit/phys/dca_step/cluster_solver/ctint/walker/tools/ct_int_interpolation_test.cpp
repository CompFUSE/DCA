// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests interpolation of G0 on the CPU.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"

#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

#include "test/unit/phys/dca_step/cluster_solver/ctint/walker/mock_parameters.hpp"

using std::cout;
using std::endl;

TEST(G0Interpolation, G0Interpolation) {
  using dca::phys::domains::time_domain;
  using dca::func::dmn_0;
  using dca::func::dmn_variadic;
  using LabelDmn = dmn_0<dca::func::dmn<2>>;

  dca::ctint::testing::MockParameters pars(M_PI, 20);

  // Initialize the domains.
  time_domain::initialize(pars);
  // dca::phys::solver::ctint::PositiveTimeDomain::initialize();

  using TestDomain = dmn_variadic<LabelDmn, dmn_0<time_domain>>;
  dca::func::function<double, TestDomain> f;

  // Write function values.
  for (int i = 0; i < time_domain::get_size(); i++) {
    const double t = time_domain::get_elements()[i];
    f(0, i) = std::sin(t);
    f(1, i) = std::sin(2 * t);
  }
  dca::phys::solver::ctint::G0Interpolation<dca::linalg::CPU> g0(f);

  for (double x : {0., 0.5, 3., M_PI - 1e-3}) {
    EXPECT_NEAR(std::sin(x), g0(x, 0), 1e-3);
    EXPECT_NEAR(std::sin(2 * x), g0(x, 1), 1e-3);

    EXPECT_NEAR(-std::sin(x), g0(-x, 0), 1e-3);
    EXPECT_NEAR(-std::sin(2 * x), -g0(x, 1), 1e-3);
  }
}
