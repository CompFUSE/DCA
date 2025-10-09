// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Author: Peter W. Doak (doakpw@ornl.gov)
//
// Type definitions for statistical tests on a square lattice.

#ifndef DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP
#define DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP

#include <string>
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/linalg/util/util_gpublas.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_mapping/cluster_exclusion.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
#include "dca/phys/dca_step/cluster_mapping/double_counting_correction.hpp"
#include "dca/phys/dca_step/cluster_mapping/update_chemical_potential.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/math/random/random.hpp"
#include "dca/phys/dca_step/cluster_mapping/update_chemical_potential.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_sp.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/testing/minimalist_printer.hpp"

namespace dca {
namespace testing {
// dca::testing::

constexpr int n_frequencies = 10;
using dca::linalg::DeviceType;
using dca::linalg::GPU;
using dca::linalg::CPU;

const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/rashba_lattice/";

const std::string default_input(test_directory + "rashba_lattice_stat_test.json");

using LatticeRashba = phys::models::RashbaHubbard<phys::domains::no_symmetry<2>>;

using dca::math::util::cutFrequency;

using dca::linalg::DeviceType;

template <typename Scalar, DeviceType DEVICE, class CONCURRENCY, class Lattice = LatticeRashba,
          ClusterSolverId solver_name = ClusterSolverId::CT_AUX, DistType DT = DistType::NONE>
struct IntegrationSetupBare {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RandomNumberGenerator = dca::math::random::StdRandomWrapper<std::mt19937_64>;
  using Concurrency = CONCURRENCY;
  using Parameters =
      phys::params::Parameters<Concurrency, Threading, profiling::NullProfiler, Model,
                               RandomNumberGenerator, solver_name,
                               dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;
  using Data = phys::DcaData<Parameters, DT>;

  /// Selector to deal with differing solver parmeter signatures
  template <ClusterSolverId name, class PARAMETERS>
  struct ClusterSolverSelector;

  template <class PARAMETERS>
  struct ClusterSolverSelector<ClusterSolverId::CT_AUX, PARAMETERS> {
    using type = dca::phys::solver::CtauxClusterSolver<DEVICE, PARAMETERS, Data>;
  };

  template <class PARAMETERS>
  struct ClusterSolverSelector<ClusterSolverId::CT_INT, PARAMETERS> {
    using type = dca::phys::solver::CtintClusterSolver<DEVICE, PARAMETERS, true>;
  };

  using QMCSolver = typename ClusterSolverSelector<solver_name, Parameters>::type;
  using ThreadedSolver = dca::phys::solver::StdThreadQmciClusterSolver<QMCSolver>;

  // Commonly used domains.
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<phys::domains::frequency_domain>;
  using LabelDomain = func::dmn_variadic<BDmn, BDmn, RDmn>;

  using KDmnDCA = dca::func::dmn_0<dca::phys::domains::cluster_domain<
      double, LatticeType::DIMENSION, dca::phys::domains::CLUSTER,
      dca::phys::domains::MOMENTUM_SPACE, dca::phys::domains::BRILLOUIN_ZONE>>;
  using KHost = dca::func::dmn_0<dca::phys::domains::cluster_domain<
      double, LatticeType::DIMENSION, dca::phys::domains::LATTICE_SP,
      dca::phys::domains::MOMENTUM_SPACE, dca::phys::domains::BRILLOUIN_ZONE>>;

  using SigmaCutDomain = dca::math::util::SigmaCutDomain<dca::math::util::details::Kdmn<>>;
  using SigmaDomain = dca::math::util::SigmaDomain<dca::math::util::details::Kdmn<>>;
  using CovarianceDomain = dca::math::util::CovarianceDomain<dca::math::util::details::Kdmn<>>;

  using coarsegraining_sp_type = dca::phys::clustermapping::CoarsegrainingSp<Parameters>;
  using lattice_map_sp_type =
      dca::phys::latticemapping::lattice_mapping_sp<Parameters, KDmnDCA, KHost>;
  using update_chemical_potential_type =
      dca::phys::clustermapping::update_chemical_potential<Parameters, Data, coarsegraining_sp_type>;

  Concurrency* concurrency_;
  Parameters parameters_;
  std::unique_ptr<Data> data_;

  std::unique_ptr<lattice_map_sp_type> lattice_mapping_obj_;
  std::unique_ptr<coarsegraining_sp_type> cluster_mapping_obj_;
  std::unique_ptr<update_chemical_potential_type> update_chemical_potential_obj_;
  IntegrationSetupBare(Concurrency* concurrency, const std::string& input_file = default_input)
      : concurrency_(concurrency), parameters_("", *concurrency_) {
    try {
      parameters_.template read_input_and_broadcast<io::JSONReader>(input_file);
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
    cluster_mapping_obj_ = std::make_unique<coarsegraining_sp_type>(parameters_);
    update_chemical_potential_obj_ =
        std::make_unique<update_chemical_potential_type>(parameters_, *data_, *cluster_mapping_obj_);
    update_chemical_potential_obj_->execute();
    data_->Sigma_cluster = data_->Sigma;
    cluster_mapping_obj_->compute_G_K_w(data_->Sigma, data_->G_k_w);
    lattice_mapping_obj_ = std::make_unique<lattice_map_sp_type>(parameters_);
    // lattice_mapping_obj_->execute(data_->Sigma, data_->Sigma_lattice_interpolated,
    //                               data_->Sigma_lattice_coarsegrained, data_->Sigma_lattice);
  }

  void performLatticeMapping() {
    lattice_mapping_obj_->execute(data_->Sigma, data_->Sigma_lattice_interpolated,
                                  data_->Sigma_lattice_coarsegrained, data_->Sigma_lattice);
  }

  void computeSigmaEtc() {
    data_.Sigma_cluster = data_.Sigma;
    data_.compute_Sigma_bands();
    data_.compute_single_particle_properties();
  }
};

}  // namespace testing
}  // namespace dca

#endif  // DCA_TEST_INTEGRATION_STATISTICAL_TESTS_SQUARE_LATTICE_SQUARE_LATTICE_SETUP_HPP
