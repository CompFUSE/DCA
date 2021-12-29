// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//

#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"

#include "gtest/gtest.h"
#include <iostream>
#include <cmath>
#include <mutex>

#include "test/unit/phys/dca_step/cluster_solver/shared_tools/interpolation/mock_parameters.hpp"

template <typename Real>
class G0InterpolationTest : public ::testing::Test {};

using FloatingPointTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(G0InterpolationTest, FloatingPointTypes);

std::once_flag flag;
TYPED_TEST(G0InterpolationTest, G0Interpolation) {
  using dca::phys::domains::time_domain;
  using dca::func::dmn_0;
  using dca::func::dmn_variadic;
  using LabelDmn = dmn_0<dca::func::dmn<2>>;

  dca::testing::MockParameters pars(M_PI, 20);

  // Initialize the domains.
  std::call_once(flag, [&] {
    time_domain::initialize(pars);
    // dca::phys::solver::PositiveTimeDomain::initialize();
  });

  using TestDomain = dmn_variadic<LabelDmn, dmn_0<time_domain>>;
  dca::func::function<double, TestDomain> f;

  // Write function values.
  for (int i = 0; i < time_domain::get_size(); i++) {
    const double t = time_domain::get_elements()[i];
    f(0, i) = std::sin(t);
    f(1, i) = std::sin(2 * t);
  }
  dca::phys::solver::G0Interpolation<dca::linalg::CPU, TypeParam> g0(f);

  for (double x : {0., 0.5, 3., M_PI - 1e-3}) {
    EXPECT_NEAR(std::sin(x), g0(x, 0), 1e-3);
    EXPECT_NEAR(std::sin(2 * x), g0(x, 1), 1e-3);

    EXPECT_NEAR(-std::sin(x), g0(-x, 0), 1e-3);
    EXPECT_NEAR(-std::sin(2 * x), -g0(x, 1), 1e-3);
  }
}
