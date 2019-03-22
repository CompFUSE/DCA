// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the coarsegraining module for single-particle functions.

#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

// This test computes the coarsegrained dispersion for a single site in a 2D square lattice model.
// Due to the periodicity of the lattice dispersion, the coarsegrained average over the entire
// Brillouin zone is zero.
TEST(CoarsegrainingSingleParticleTest, ComputeCoarsegrainedDispersion) {
  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;
  using Model = phys::models::TightBindingModel<Lattice>;

  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, parallel::stdthread, profiling::NullProfiler, Model,
                               void /*RandomNumberGenerator*/, phys::solver::CT_AUX>;

  using CoarsegrainingType = phys::clustermapping::CoarsegrainingSp<ParametersType>;

  using KClusterDmn = CoarsegrainingType::KClusterDmn;
  using NuDmn = CoarsegrainingType::NuDmn;
  using ComplexType = CoarsegrainingType::Complex;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);

  parameters.set_t(1.);
  parameters.set_k_mesh_recursion(3);
  parameters.set_coarsegraining_periods(0);
  parameters.set_quadrature_rule(3);
  parameters.set_cluster({{0, 1}, {1, 0}});

  parameters.update_model();
  parameters.update_domains();

  CoarsegrainingType coarsegrain_sp(parameters);

  func::function<ComplexType, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>> eps_K_cg;

  coarsegrain_sp.computeCoarsegrainedDispersion(eps_K_cg);

  // The dispersion is purely real.
  for (int i = 0; i < eps_K_cg.size(); ++i) {
    EXPECT_DOUBLE_EQ(0., eps_K_cg(i).imag());
  }

  // Off-diagonal elements should be zero.
  EXPECT_DOUBLE_EQ(0., eps_K_cg(0, 1, 0).real());
  EXPECT_DOUBLE_EQ(0., eps_K_cg(1, 0, 0).real());

  // The coarsegrained average over the entire Brillouin zone shuld be zero.
  EXPECT_LE(eps_K_cg(0, 0, 0).real(), 1.e-10);
  EXPECT_LE(eps_K_cg(1, 1, 0).real(), 1.e-10);
}
