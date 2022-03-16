// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// Type definitions for statistical tests on NiO material model

#ifndef DCA_TEST_INTEGRATION_STATISTICAL_TESTS_NIO_SETUP_HPP
#define DCA_TEST_INTEGRATION_STATISTICAL_TESTS_NIO_SETUP_HPP

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp"



#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/math/random/random.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"

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

#ifdef DCA_HAVE_GPU
constexpr DeviceType default_device = GPU;
#else
constexpr DeviceType default_device = CPU;
#endif  // DCA_HAVE_CUDA

const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/real_materials/";

using Model = dca::phys::models::TightBindingModel<dca::phys::models::material_lattice<
      dca::phys::models::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>>;
using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::ranlux48_base>;

using dca::ClusterSolverId;

template <ClusterSolverId name>
using TestParameters =
    dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                  Threading, dca::profiling::NullProfiler, Model,
                                  RandomNumberGenerator, name>;

template <ClusterSolverId name>
using DcaData = dca::phys::DcaData<TestParameters<name>>;

template <ClusterSolverId name, DeviceType device>
struct ClusterSolverSelector;

template <DeviceType device>
struct ClusterSolverSelector<dca::ClusterSolverId::SS_CT_HYB, device> {
  using type = dca::phys::solver::SsCtHybClusterSolver<device, TestParameters<dca::ClusterSolverId::SS_CT_HYB>, DcaData<dca::ClusterSolverId::SS_CT_HYB>>;
};
template <ClusterSolverId name, DeviceType device>
using QuantumClusterSolver = typename ClusterSolverSelector<name, device>::type;

template <ClusterSolverId name, DeviceType device>
using ThreadedSolver =
    dca::phys::solver::StdThreadQmciClusterSolver<QuantumClusterSolver<name, device>>;

using dca::func::dmn_0;
using dca::func::dmn_variadic;

// using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn>;
// using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn>;
// using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn>;
using dca::math::util::cutFrequency;

}  // namespace testing
}  // namespace dca

#endif  // DCA_TEST_INTEGRATION_STATISTICAL_TESTS_NIO_SETUP_HPP
