// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Type definitions for statistical tests on a square lattice.

#ifndef DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP
#define DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/math/random/random.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"

namespace dca {
namespace testing {
// dca::testing::

const int n_frequencies = 10;

const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/square_lattice/";

using Model =
    dca::phys::models::TightBindingModel<dca::phys::models::square_lattice<dca::phys::domains::D4>>;
using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using ParametersType =
    dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                  Threading, dca::profiling::NullProfiler, Model,
                                  RandomNumberGenerator, dca::phys::solver::CT_AUX>;
using DcaData = dca::phys::DcaData<ParametersType>;
using QuantumClusterSolver =
    dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, ParametersType, DcaData>;
using ThreadedSolver = dca::phys::solver::StdThreadQmciClusterSolver<QuantumClusterSolver>;

using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn>;
using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn>;
using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn>;
using dca::math::util::cutFrequency;

}  // testing
}  // dca

#endif  // DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP
