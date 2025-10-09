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

#include "dca/config/threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"

#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/math/random/random.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"

#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

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

using Scalar = double;

const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/real_materials/";

using Model = dca::phys::models::TightBindingModel<dca::phys::models::material_lattice<
  dca::phys::models::Material::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>>;
using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::ranlux48_base>;

using dca::ClusterSolverId;

template <class Concurrency, ClusterSolverId name>
using TestParameters =
    dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler, Model,
                                  RandomNumberGenerator, name, dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;

template <ClusterSolverId name, class Concurrency>
using DcaData = dca::phys::DcaData<TestParameters<Concurrency, name>>;

template <ClusterSolverId name, class Concurrency, DeviceType device>
struct ClusterSolverSelector;

template <class Concurrency, DeviceType device>
struct ClusterSolverSelector<ClusterSolverId::CT_INT, Concurrency, device> {
  using type = dca::phys::solver::CtintClusterSolver<
      device, TestParameters<Concurrency, ClusterSolverId::CT_INT>, true, DistType::NONE>;
};

template <class Concurrency, DeviceType device>
struct ClusterSolverSelector<dca::ClusterSolverId::SS_CT_HYB, Concurrency, device> {
  using type = dca::phys::solver::SsCtHybClusterSolver<
      device, TestParameters<Concurrency, dca::ClusterSolverId::SS_CT_HYB>,
    DcaData<dca::ClusterSolverId::SS_CT_HYB, Concurrency>, DistType::NONE>;
};

template <ClusterSolverId name, class Concurrency, DeviceType device>
using QuantumClusterSolver = typename ClusterSolverSelector<name, Concurrency, device>::type;

template <ClusterSolverId name, class Concurrency, DeviceType device>
using ThreadedSolver =
    dca::phys::solver::StdThreadQmciClusterSolver<QuantumClusterSolver<name, Concurrency, device>>;

using dca::func::dmn_0;
using dca::func::dmn_variadic;

// using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn>;
// using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn>;
// using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn>;
using dca::math::util::cutFrequency;

}  // namespace testing
}  // namespace dca

#endif  // DCA_TEST_INTEGRATION_STATISTICAL_TESTS_NIO_SETUP_HPP
