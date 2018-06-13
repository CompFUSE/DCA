// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests wannier_interpolation.hpp.

#include "dca/phys/domains/cluster/interpolation/wannier_interpolation/wannier_interpolation.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

TEST(WannierInterpolationTest, TightBindingHamiltonianSquareLattice) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;
  using Model = phys::models::TightBindingModel<Lattice>;

  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, parallel::NoThreading, profiling::NullProfiler,
                               Model, void /*RandomNumberGenerator*/, phys::solver::CT_AUX>;

  using b = func::dmn_0<phys::domains::electron_band_domain>;
  using s = func::dmn_0<phys::domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;

  // Defined by input parameter "cluster".
  using k_DCA = func::dmn_0<
      phys::domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, phys::domains::CLUSTER,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;
  // Defined by input parameter "H(k) grid-size".
  using k_LDA = func::dmn_0<
      phys::domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, phys::domains::LATTICE_SP,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::PARALLELLEPIPEDUM>>;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);
  parameters.read_input_and_broadcast<io::JSONReader>(
      DCA_SOURCE_DIR
      "/test/unit/phys/domains/cluster/interpolation/wannier_interpolation/input.json");
  parameters.update_model();
  parameters.update_domains();

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA>> H_DCA_interpolated;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_DCA>> H_DCA_direct;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_LDA>> H_LDA;

  phys::models::square_lattice<PointGroup>::initialize_H_0(parameters, H_LDA);

  // Compute H_DCA by interpolating H_LDA.
  phys::domains::wannier_interpolation<k_LDA, k_DCA>::execute(H_LDA, H_DCA_interpolated);

  // Directly compute H_DCA.
  phys::models::square_lattice<PointGroup>::initialize_H_0(parameters, H_DCA_direct);

  for (int b1 = 0; b1 < b::dmn_size(); ++b1)
    for (int s1 = 0; s1 < s::dmn_size(); ++s1)
      for (int b2 = 0; b2 < b::dmn_size(); ++b2)
        for (int s2 = 0; s2 < s::dmn_size(); ++s2)
          for (int k = 0; k < k_DCA::dmn_size(); ++k) {
            EXPECT_NEAR(H_DCA_interpolated(b1, s1, b2, s2, k).real(),
                        H_DCA_direct(b1, s1, b2, s2, k).real(), 1.e-13);
            EXPECT_NEAR(H_DCA_interpolated(b1, s1, b2, s2, k).imag(),
                        H_DCA_direct(b1, s1, b2, s2, k).imag(), 1.e-13);
          }
}
