// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This tests space_transform_2D.hpp

#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include <string>
#include <random>

#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#ifdef DCA_HAVE_ADIOS2
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#endif
dca::parallel::NoConcurrency* concurrency_ptr;

using Model = dca::phys::models::TightBindingModel<
    dca::phys::models::twoband_chain<dca::phys::domains::no_symmetry<2>>>;
template <typename SCALAR>
using NumTraits = dca::NumericalTraits<dca::util::RealAlias<SCALAR>, SCALAR>;

using Concurrency = dca::parallel::NoConcurrency;

const std::string input_dir = DCA_SOURCE_DIR "/test/unit/math/function_transform/";

using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using SDmn = dca::func::dmn_0<dca::phys::domains::electron_spin_domain>;
using WPosDmn =
    dca::func::dmn_0<dca::phys::domains::vertex_frequency_domain<dca::phys::domains::COMPACT_POSITIVE>>;
using WDmn =
    dca::func::dmn_0<dca::phys::domains::vertex_frequency_domain<dca::phys::domains::COMPACT>>;

template <typename SCALAR>
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  Model, void, dca::ClusterSolverId::CT_AUX, NumTraits<SCALAR>>;

// Perform the test in double and single precision.
template <typename T>
class SpaceTransform2DTest : public ::testing::Test {};
using TestTypes = ::testing::Types<float, double>;  //, std::complex<double>>;
TYPED_TEST_CASE(SpaceTransform2DTest, TestTypes);

TYPED_TEST(SpaceTransform2DTest, Execute) {
  using Scalar = TypeParam;

  using KDmn = typename Parameters<Scalar>::KClusterDmn;
  using RDmn = typename Parameters<Scalar>::RClusterDmn;

  Parameters<Scalar> pars("", *concurrency_ptr);
  pars.template read_input_and_broadcast<dca::io::JSONReader>(input_dir + "input.json");
  pars.update_model();
  pars.update_domains();

  using dca::func::dmn_variadic;
  using dca::func::function;
  using Real = dca::util::RealAlias<Scalar>;
  using Complex = std::complex<Real>;
  std::mt19937_64 rng(0);
  std::uniform_real_distribution<Real> distro(-1, 1);

  // Initialize the input function randomly.
  function<Complex, dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>> f_in;
  for (int i = 0; i < f_in.size(); ++i)
    f_in(i) = Complex(distro(rng), distro(rng));
  auto f_in_cpy = f_in;

  dca::func::function<Complex, dca::func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, WPosDmn, WDmn>> f_out;
  dca::math::transform::SpaceTransform2D<RDmn, KDmn, BDmn, SDmn, Scalar>::execute(f_in_cpy, f_out);

  const auto im = std::complex<double>(0, 1);

  const auto& k = KDmn::get_elements();
  const auto& r = RDmn::get_elements();

  // Displacement vectors.
  std::vector<std::vector<double>> a;
  for (const auto& el : BDmn::get_elements())
    a.push_back(el.a_vec);

  const int nb = BDmn::dmn_size();
  const int nc = RDmn::dmn_size();
  const int nw = WPosDmn::dmn_size();

  for (int w2 = 0; w2 < 2 * nw; ++w2)
    for (int w1 = 0; w1 < nw; ++w1)
      for (int k2 = 0; k2 < nc; ++k2)
        for (int k1 = 0; k1 < nc; ++k1)
          for (int s = 0; s < 2; ++s)
            for (int b2 = 0; b2 < nb; ++b2)
              for (int b1 = 0; b1 < nb; ++b1) {
                // Compute the expected value.
                Complex expected = 0;
                using namespace dca::math::util;  // vector operations.
                for (int r2 = 0; r2 < nc; ++r2)
                  for (int r1 = 0; r1 < nc; ++r1)
                    expected += static_cast<Complex>(std::exp(
                                    im * (k[k1] * (r[r1] + a[b1]) - k[k2] * (r[r2] + a[b2])))) *
                                f_in(r1, r2, b1, b2, s, w1, w2);
                expected /= nc;

                EXPECT_NEAR(expected.real(), f_out(b1, b2, s, k1, k2, w1, w2).real(), 1e-5);
                EXPECT_NEAR(expected.imag(), f_out(b1, b2, s, k1, k2, w1, w2).imag(), 1e-5);
              }
}

int main(int argc, char** argv) {
#ifdef DCA_HAVE_ADIOS2
  // Until I figure out how to prevent globals inside the adios library from calling
  // MPI_COMM_dup why're stuck with this if it is linked.
  dca::parallel::MPIConcurrency mpi_concurrency(argc, argv);
#endif
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;

  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
