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

#ifndef DCA_TEST_INTEGRATION_STATISTICAL_TESTS_FEAS_FE_AS_SETUP_HPP
#define DCA_TEST_INTEGRATION_STATISTICAL_TESTS_FEAS_FE_AS_SETUP_HPP

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"
#include "dca/math/random/random.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_Cu.hpp"
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

constexpr int n_frequencies = 10;
using dca::linalg::DeviceType;
using dca::linalg::GPU;
using dca::linalg::CPU;

#ifdef DCA_HAVE_CUDA
constexpr DeviceType default_device = GPU;
#else
constexpr DeviceType default_device = CPU;
#endif  // DCA_HAVE_CUDA

const std::string test_directory = DCA_SOURCE_DIR "/test/integration/statistical_tests/FeAs/";

using Model =
    dca::phys::models::TightBindingModel<dca::phys::models::FeAsLattice<dca::phys::domains::D4>>;
using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::mt19937_64>;

using dca::phys::solver::CT_INT;

using ParametersType =
    dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                  dca::parallel::stdthread, dca::profiling::NullProfiler, Model,
                                  RandomNumberGenerator, CT_INT>;

using DcaData = dca::phys::DcaData<ParametersType>;

template <DeviceType device>
using QuantumClusterSolver = dca::phys::solver::CtintClusterSolver<device, ParametersType, true>;

template <DeviceType device>
using ThreadedSolver = dca::phys::solver::StdThreadQmciClusterSolver<QuantumClusterSolver<device>>;

using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn>;
using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn>;
using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn>;
using dca::math::util::cutFrequency;

}  // namespace testing
}  // namespace dca

#endif  // DCA_TEST_INTEGRATION_STATISTICAL_TESTS_FEAS_FE_AS_SETUP_HPP
