
#include "dca/config/haves_defines.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation_gpu.hpp"

#include <cmath>
#include <mutex>
#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/interpolation/mock_parameters.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/kernels_interface.hpp"
#include "dca/linalg/util/gpu_type_mapping.hpp"
#include "dca/util/type_utils.hpp"
using std::cout;
using std::endl;

template <typename ScalarType>
class G0InterpolationGpuTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double, std::complex<float>, std::complex<double>>;
TYPED_TEST_CASE(G0InterpolationGpuTest, TestTypes);

using dca::util::RealAlias;
using dca::util::IsComplex_t;

std::once_flag flag;
TYPED_TEST(G0InterpolationGpuTest, G0Interpolation) {
  using Scalar = TypeParam;
  using Real = RealAlias<Scalar>;
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
  dca::func::function<Scalar, TestDomain> f;

  // Write function values.
  for (int i = 0; i < time_domain::get_size(); i++) {
    const Real t = time_domain::get_elements()[i];
    if constexpr(IsComplex_t<Scalar>::value) {
      f(0, i) = {std::exp(t), -std::sin(t)};
      f(1, i) = {std::sin(t), -std::exp(t)};
    }
    else {
    f(0, i) = std::exp(t);
    f(1, i) = std::sin(t);
    }
  }
  dca::phys::solver::G0Interpolation<dca::linalg::CPU, Scalar> g0_cpu(f);

  dca::phys::solver::G0Interpolation<dca::linalg::GPU, Scalar> g0_gpu(f);

  constexpr Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
  for (Real x : {0., 0.5, 3., M_PI - 1e-3}) {
    EXPECT_NEAR(std::real(g0_cpu(x, 0)), std::real(g0_gpu(x, 0)), tolerance);
    EXPECT_NEAR(std::real(g0_cpu(x, 1)), std::real(g0_gpu(x, 1)), tolerance);
    if constexpr(IsComplex_t<Scalar>::value) {
      EXPECT_NEAR(std::imag(g0_cpu(x, 0)), std::imag(g0_gpu(x, 0)), tolerance);
      EXPECT_NEAR(std::imag(g0_cpu(x, 1)), std::imag(g0_gpu(x, 1)), tolerance);
    }
  }
}
