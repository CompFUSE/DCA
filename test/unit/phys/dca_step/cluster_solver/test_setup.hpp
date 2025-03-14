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
#include <iostream>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/json/json_reader.hpp"

#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/interaction_vertices.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"
#include "dca/phys/models/analytic_hamiltonians/Moire_Hubbard.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/models/analytic_hamiltonians/Kagome_hubbard.hpp"
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
using LatticeHund = phys::models::HundLattice<phys::domains::D4>;
using LatticeKagome = phys::models::KagomeHubbard<phys::domains::no_symmetry<2>>;
using LatticeRashba = phys::models::RashbaHubbard<phys::domains::no_symmetry<2>>;
using LatticeMoireHubbard = phys::models::moire_hubbard<phys::domains::no_symmetry<2>>;

template <typename Scalar, class Lattice = LatticeSquare,
          ClusterSolverId solver_name = ClusterSolverId::CT_AUX,
          const char* input_name = default_input, DistType DT = DistType::NONE>
struct G0Setup : public ::testing::Test {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RngType = testing::StubRng;
  using Concurrency = parallel::NoConcurrency;
  using Parameters =
      phys::params::Parameters<Concurrency, parallel::NoThreading, profiling::NullProfiler, Model, RngType,
                               solver_name, dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
  using Data = phys::DcaData<Parameters, DT>;

  // Commonly used domains.
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<phys::domains::frequency_domain>;
  using LabelDomain = func::dmn_variadic<BDmn, BDmn, RDmn>;

  Concurrency concurrency_;
  Parameters parameters_;
  std::unique_ptr<Data> data_;

  G0Setup() : concurrency_(0, nullptr), parameters_("", concurrency_) {}

  virtual void SetUp() {
    try {
      parameters_.template read_input_and_broadcast<io::JSONReader>(input_name);
    }
    catch (const std::exception& r_w) {
      throw std::runtime_error(r_w.what());
    }
    catch (...) {
      throw std::runtime_error("Input parsing failed!");
    }
    parameters_.update_model();
    static bool domain_initialized = false;
    if (!domain_initialized) {
      parameters_.update_domains();
      domain_initialized = true;
    }
    data_ = std::make_unique<Data>(parameters_);
    data_->initialize();
  }

  virtual void TearDown() {}

  auto& getParameters() {
      return parameters_;
  }

};

template <typename Scalar, class Lattice = LatticeSquare,
          ClusterSolverId solver_name = ClusterSolverId::CT_AUX,
          const char* input_name = default_input, DistType DT = DistType::NONE>
struct G0SetupBare {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RngType = testing::StubRng;
  using Concurrency = parallel::NoConcurrency;
  using Parameters =
      phys::params::Parameters<Concurrency, parallel::NoThreading, profiling::NullProfiler, Model, RngType,
                               solver_name, dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
  using Data = phys::DcaData<Parameters, DT>;

  // Commonly used domains.
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<phys::domains::frequency_domain>;
  using LabelDomain = func::dmn_variadic<BDmn, BDmn, RDmn>;

  Concurrency concurrency_;
  Parameters parameters_;
  std::unique_ptr<Data> data_;

  G0SetupBare() : concurrency_(0, nullptr), parameters_("", concurrency_) {}

  void SetUp() {
    try {
      parameters_.template read_input_and_broadcast<io::JSONReader>(input_name);
    }
    catch (const std::exception& r_w) {
      throw std::runtime_error(r_w.what());
    }
    catch (...) {
      throw std::runtime_error("Input parsing failed!");
    }
    parameters_.update_model();

    parameters_.update_domains();

    data_ = std::make_unique<Data>(parameters_);
    data_->initialize();
  }

  void TearDown() {}

};

}  // namespace testing
}  // namespace dca

#endif  // DCA_TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_TEST_SETUP_HPP
