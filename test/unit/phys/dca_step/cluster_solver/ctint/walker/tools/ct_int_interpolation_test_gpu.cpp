#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation_gpu.hpp"

#include <cmath>
#include <mutex>
#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/ctint/device_helper/ctint_helper.cuh"
#include "test/unit/phys/dca_step/cluster_solver/ctint/walker/mock_parameters.hpp"

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

  dca::ctint::testing::MockParameters pars(M_PI, 20);

  // Initialize the domains.
  std::call_once(flag, [&] {
    time_domain::initialize(pars);
    // dca::phys::solver::ctint::PositiveTimeDomain::initialize();
  });

  using TestDomain = dmn_variadic<LabelDmn, dmn_0<time_domain>>;
  dca::func::function<double, TestDomain> f;

  // Write function values.
  for (int i = 0; i < time_domain::get_size(); i++) {
    const Real t = time_domain::get_elements()[i];
    f(0, i) = std::exp(t);
    f(1, i) = std::sin(t);
  }
  dca::phys::solver::ctint::G0Interpolation<dca::linalg::CPU, Real> g0_cpu(f);

  dca::phys::solver::ctint::G0Interpolation<dca::linalg::GPU, Real> g0_gpu(f);

  constexpr Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
  for (Real x : {0., 0.5, 3., M_PI - 1e-3}) {
    EXPECT_NEAR(g0_cpu(x, 0), g0_gpu(x, 0), tolerance);
    EXPECT_NEAR(g0_cpu(x, 1), g0_gpu(x, 1), tolerance);
  }
}
