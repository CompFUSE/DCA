// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the FunctionTransform class for the real to momentum space transformations.

#include "dca/math/function_transform/function_transform.hpp"

#include "gtest/gtest.h"
#include <string>

#include "dca/io/json/json_reader.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"

#ifdef DCA_HAVE_ADIOS2
adios2::ADIOS* adios_ptr;
#endif
#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
dca::parallel::MPIConcurrency* concurrency_ptr;
#else
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
dca::parallel::NoConcurrency* concurrency_ptr;
#endif

using Model =
    dca::phys::models::TightBindingModel<dca::phys::models::bilayer_lattice<dca::phys::domains::D4>>;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  Model, void, dca::ClusterSolverId::CT_AUX>;

const std::string input_dir = DCA_SOURCE_DIR "/test/unit/math/function_transform/";

using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
using SDmn = dca::func::dmn_0<dca::phys::domains::electron_spin_domain>;
using NuDmn = dca::func::dmn_variadic<BDmn, SDmn>;
using WDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
using KDmn = Parameters::KClusterDmn;
using RDmn = Parameters::RClusterDmn;

const std::vector<std::vector<double>> a_vecs{std::vector<double>{0, 0},
                                              std::vector<double>{0.25, 0.25}};

class FunctionTransformTest : public ::testing::Test {};

void initialize() {
  static bool initialized = false;
  if (!initialized) {
    Concurrency concurrency(0, nullptr);
    Parameters pars("", concurrency);
    pars.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "input.json");
    pars.update_model();
    pars.update_domains();
    BDmn::parameter_type::setAVectors(a_vecs);
    initialized = true;
  }
}

template <class InpDmn, class OutDmn>
void spTestImplementation(const bool direct) {
  using namespace dca::func;
  using Real = double;
  using Complex = std::complex<Real>;
  function<Complex, dmn_variadic<NuDmn, NuDmn, InpDmn, WDmn>> f_b_b_in;
  function<Complex, dmn_variadic<InpDmn, WDmn>> f_in;

  // Initialize the input function.
  for (int w = 0; w < Parameters::WDmn::dmn_size(); ++w)
    for (int r = 0; r < InpDmn::dmn_size(); ++r) {
      const Complex val(r * r + 0.5 * r, w);
      f_in(r, w) = val;
      for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
        for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1)
          f_b_b_in(b1, 0, b2, 0, r, w) = f_b_b_in(b1, 1, b2, 1, r, w) = val;
    }

  function<Complex, dmn_variadic<NuDmn, NuDmn, OutDmn, WDmn>> f_b_b_out;
  function<Complex, dmn_variadic<OutDmn, WDmn>> f_out;

  // Regular transform.
  dca::math::transform::FunctionTransform<InpDmn, OutDmn>::execute(f_in, f_out);
  // Transform with phase factor.
  dca::math::transform::FunctionTransform<InpDmn, OutDmn>::execute(f_b_b_in, f_b_b_out);

  const int exp_sign = direct ? 1 : -1;
  const double norm = direct ? 1 : 1. / InpDmn::dmn_size();

  constexpr double tolerance = 1e-10;

  for (int w = 0; w < WDmn::dmn_size(); ++w) {
    // Test regular transform.
    for (int o = 0; o < OutDmn::dmn_size(); ++o) {
      const auto& o_dmn_val = OutDmn::get_elements()[o];
      Complex val(0);

      for (int i = 0; i < InpDmn::dmn_size(); ++i) {
        const auto& i_dmn_val = InpDmn::get_elements()[i];
        val += f_in(i, w) *
               std::exp(Complex(0., exp_sign * dca::math::util::innerProduct(i_dmn_val, o_dmn_val)));
      }
      val *= norm;

      EXPECT_LE(std::abs(val - f_out(o, w)), tolerance);

      // Test transform with phase factors.
      for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1) {
        for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2) {
          const auto a_diff = dca::math::util::subtract(a_vecs[b2], a_vecs[b1]);
          for (int o = 0; o < OutDmn::dmn_size(); ++o) {
            const auto& o_dmn_val = OutDmn::get_elements()[o];
            Complex val(0);

            for (int i = 0; i < InpDmn::dmn_size(); ++i) {
              const auto& i_dmn_val = InpDmn::get_elements()[i];
              const auto& k_val = direct ? o_dmn_val : i_dmn_val;
              const Complex phase_factor =
                  std::exp(Complex(0, exp_sign * dca::math::util::innerProduct(k_val, a_diff)));
              val += f_b_b_in(b1, 0, b2, 0, i, w) * phase_factor *
                     std::exp(Complex(
                         0., exp_sign * dca::math::util::innerProduct(i_dmn_val, o_dmn_val)));
            }
            val *= norm;

            EXPECT_LE(std::abs(val - f_b_b_out(b1, 0, b2, 0, o, w)), tolerance);
          }
        }
      }
    }
  }
}

TEST(FunctionTransformTest, SpaceToMomentumCmplx) {
  spTestImplementation<RDmn, KDmn>(true);
}

TEST(FunctionTransformTest, MomentumToSpaceCmplx) {
  spTestImplementation<KDmn, RDmn>(false);
}

int main(int argc, char** argv) {
#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#else
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#endif

#ifdef DCA_HAVE_ADIOS2
  //ADIOS expects MPI_COMM pointer or nullptr
  adios2::ADIOS adios("", concurrency_ptr->get(), "C++");
  adios_ptr = &adios;
#endif

  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  initialize();

  int result = RUN_ALL_TESTS();
  return result;
}
