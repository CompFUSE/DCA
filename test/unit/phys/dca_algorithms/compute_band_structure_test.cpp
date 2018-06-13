// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests compute_band_structure.hpp.

#include "dca/phys/dca_algorithms/compute_band_structure.hpp"

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/brillouin_zone_cut_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

TEST(ComputeBandStructureTest, Execute) {
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
  using k_domain_cut_dmn_type = func::dmn_0<phys::domains::brillouin_zone_cut_domain<101>>;
  using nu_k_cut = func::dmn_variadic<nu, k_domain_cut_dmn_type>;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);
  parameters.read_input_and_broadcast<io::JSONReader>(DCA_SOURCE_DIR
                                                      "/test/unit/phys/dca_algorithms/input.json");
  parameters.update_model();
  parameters.update_domains();

  func::function<double, nu_k_cut> band_structure;
  phys::compute_band_structure::execute(parameters, band_structure);

  // Check spin symmetry.
  for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
    for (int k_ind = 0; k_ind < k_domain_cut_dmn_type::dmn_size(); ++k_ind)
      EXPECT_DOUBLE_EQ(band_structure(b_ind, 0, k_ind), band_structure(b_ind, 1, k_ind));

  // Check min and max elements for one spin type.
  std::vector<double> spin_up_values;
  for (int k_ind = 0; k_ind < k_domain_cut_dmn_type::dmn_size(); ++k_ind)
    spin_up_values.push_back(band_structure(0, 0, k_ind));

  EXPECT_DOUBLE_EQ(-6., *(std::min_element(spin_up_values.begin(), spin_up_values.end())));
  EXPECT_DOUBLE_EQ(2., *(std::max_element(spin_up_values.begin(), spin_up_values.end())));

  // Check some values.
  EXPECT_DOUBLE_EQ(-6., band_structure(0, 0, 0));
  EXPECT_DOUBLE_EQ(-5.9961307263015469, band_structure(0, 0, 1));
  EXPECT_DOUBLE_EQ(-5.9845322608976534, band_structure(0, 0, 2));
  EXPECT_DOUBLE_EQ(-5.6207103439591695, band_structure(0, 0, 10));
  EXPECT_DOUBLE_EQ(-0.062690965389416001, band_structure(0, 0, 50));
  EXPECT_DOUBLE_EQ(1.999999532034358, band_structure(0, 0, 100));
  EXPECT_DOUBLE_EQ(2., band_structure(0, 0, 200));
  EXPECT_DOUBLE_EQ(1.9826352482222638, band_structure(0, 0, 300));
  EXPECT_DOUBLE_EQ(-5.9690794894531063, band_structure(0, 0, 400));
  EXPECT_DOUBLE_EQ(-5.9980651291679523, band_structure(0, 0, 403));
}
