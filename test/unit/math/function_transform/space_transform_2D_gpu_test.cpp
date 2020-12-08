// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This test confronts a 2D space to momentum function transform executed on the GPU, with the same
// transform executed on the CPU.

#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"

#include "gtest/gtest.h"
#include <string>

#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"

using Model = dca::phys::models::TightBindingModel<
    dca::phys::models::twoband_chain<dca::phys::domains::no_symmetry<2>>>;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  Model, void, dca::phys::solver::CT_AUX>;

const std::string input_dir = DCA_SOURCE_DIR "/test/unit/math/function_transform/";

using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using SDmn = dca::func::dmn_0<dca::phys::domains::electron_spin_domain>;
using KDmn = typename Parameters::KClusterDmn;
using RDmn = typename Parameters::RClusterDmn;
using WPosDmn =
    dca::func::dmn_0<dca::phys::domains::vertex_frequency_domain<dca::phys::domains::COMPACT_POSITIVE>>;
using WDmn =
    dca::func::dmn_0<dca::phys::domains::vertex_frequency_domain<dca::phys::domains::COMPACT>>;

void initialize() {
  static bool initialized = false;
  if (!initialized) {
    Concurrency concurrency(0, nullptr);
    Parameters pars("", concurrency);
    pars.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "input.json");
    pars.update_model();
    pars.update_domains();

    dca::linalg::util::initializeMagma();
    initialized = true;
  }
}

template <typename Scalar, dca::linalg::DeviceType device>
using RMatrix = dca::linalg::ReshapableMatrix<Scalar, device>;

template <typename Real>
using SpaceTransform2DGpuTest = ::testing::Test;
using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(SpaceTransform2DGpuTest, TestTypes);

TYPED_TEST(SpaceTransform2DGpuTest, Execute) {
  initialize();

  using dca::func::dmn_variadic;
  using dca::func::function;
  using Real = TypeParam;
  using Complex = std::complex<Real>;
  function<Complex, dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>> f_in;
  dca::linalg::ReshapableMatrix<Complex, dca::linalg::CPU> M_in;

  // Initialize the input function.
  const int nb = BDmn::dmn_size();
  const int nr = RDmn::dmn_size();
  const int nw = WPosDmn::dmn_size();
  M_in.resizeNoCopy(std::make_pair(nb * nr * nw, nb * nr * 2 * nw));
  for (int w2 = 0; w2 < 2 * nw; ++w2)
    for (int w1 = 0; w1 < nw; ++w1)
      for (int r2 = 0; r2 < nr; ++r2)
        for (int r1 = 0; r1 < nr; ++r1)
          for (int b2 = 0; b2 < nb; ++b2)
            for (int b1 = 0; b1 < nb; ++b1) {
              // Initialize the input with some arbitrary function of the indices.
              const Complex val(r1 * r1 + b1 - 0.5 * w1, r2 * r2 + b2 - 0.5 * w2);

              auto index = [=](int r, int b, int w) { return r + nr * b + nb * nr * w; };
              f_in(r1, r2, b1, b2, 0, w1, w2) = M_in(index(r1, b1, w1), index(r2, b2, w2)) = val;
            }

  // Transform on the CPU.
  function<Complex, dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, WPosDmn, WDmn>> f_out;
  dca::math::transform::SpaceTransform2D<RDmn, KDmn, Real>::execute(f_in, f_out);

  // Transform on the GPU.
  dca::linalg::ReshapableMatrix<Complex, dca::linalg::GPU, dca::config::McOptions::TpAllocator<Complex>>
      M_dev(M_in);

  dca::linalg::util::MagmaQueue queue;

  dca::math::transform::SpaceTransform2DGpu<RDmn, KDmn, Real> transform_obj(nw, queue);
  transform_obj.execute(M_dev);

  queue.getStream().sync();

  constexpr Real tolerance = std::numeric_limits<Real>::epsilon() * 500;

  RMatrix<Complex, dca::linalg::CPU> M_out(M_dev);
  for (int w2 = 0; w2 < 2 * nw; ++w2)
    for (int w1 = 0; w1 < nw; ++w1)
      for (int r2 = 0; r2 < nr; ++r2)
        for (int r1 = 0; r1 < nr; ++r1)
          for (int b2 = 0; b2 < nb; ++b2)
            for (int b1 = 0; b1 < nb; ++b1) {
              const Complex val1 = f_out(b1, b2, 0, r1, r2, w1, w2);
              auto index = [=](int b, int r, int w) { return b + nb * r + nb * nr * w; };
              const Complex val2 = M_out(index(b1, r1, w1), index(b2, r2, w2));

              EXPECT_LE(std::abs(val1 - val2), tolerance);
            }
}
