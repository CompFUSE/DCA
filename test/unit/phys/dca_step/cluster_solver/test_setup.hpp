// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a setup for Parameters and DcaData used by cluster solver tests.

#ifndef DCA_TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_TEST_SETUP_HPP
#define DCA_TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_TEST_SETUP_HPP

#include <memory>
#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"

namespace dca {
namespace testing {
// dca::testing::

constexpr char default_input[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/input.json";

using LatticeSquare = phys::models::square_lattice<phys::domains::D4>;
using LatticeBilayer = phys::models::bilayer_lattice<phys::domains::D4>;

template <class Lattice = LatticeSquare, phys::solver::ClusterSolverName solver_name = phys::solver::CT_AUX,
          const char* input_name = default_input>
struct G0Setup : public ::testing::Test {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RngType = testing::StubRng;
  using Concurrency = parallel::NoConcurrency;
  using Parameters = phys::params::Parameters<Concurrency, parallel::NoThreading,
                                              profiling::NullProfiler, Model, RngType, solver_name>;
  using Data = phys::DcaData<Parameters>;

  // Commonly used domains.
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<phys::domains::frequency_domain>;

  Concurrency concurrency;
  Parameters parameters_;
  std::unique_ptr<Data> data_;

  G0Setup() : concurrency(0, nullptr), parameters_("", concurrency) {}

  virtual void SetUp() {
    parameters_.template read_input_and_broadcast<io::JSONReader>(input_name);

    parameters_.update_model();
    static bool domain_initialized = false;
    if (!domain_initialized) {
      parameters_.update_domains();
      domain_initialized = true;
    }
    data_= std::make_unique<Data>(parameters_);
    data_->initialize();
  }

  virtual void TearDown() {}
};

}  // testing
}  // dca

#endif  // DCA_TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_TEST_SETUP_HPP
