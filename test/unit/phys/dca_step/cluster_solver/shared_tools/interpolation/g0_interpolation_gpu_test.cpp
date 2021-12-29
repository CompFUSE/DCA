// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//

#include <cmath>
#include <mutex>
#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/interpolation/mock_parameters.hpp"

using std::cout;
using std::endl;

template <typename ScalarType>
class G0InterpolationGpuTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(G0InterpolationGpuTest, TestTypes);

std::once_flag flag;
TYPED_TEST(G0InterpolationGpuTest, G0Interpolation) {
  using Real = TypeParam;

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
    const Real t = time_domain::get_elements()[i];
    f(0, i) = std::exp(t);
    f(1, i) = std::sin(t);
  }
  dca::phys::solver::G0Interpolation<dca::linalg::CPU, Real> g0_cpu(f);

  dca::phys::solver::G0Interpolation<dca::linalg::GPU, Real> g0_gpu(f);

  constexpr Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
  for (Real x : {0., 0.5, 3., M_PI - 1e-3}) {
    EXPECT_NEAR(g0_cpu(x, 0), g0_gpu(x, 0), tolerance);
    EXPECT_NEAR(g0_cpu(x, 1), g0_gpu(x, 1), tolerance);
  }
}
