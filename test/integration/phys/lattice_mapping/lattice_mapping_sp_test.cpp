// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test of / example for the single-particle (sp) lattice mapping.
// The value-parameterized test is run with
// - a DCA+ cluster self-energy for which a lattice self-energy can be found with the given
//   tolerance and maximum number of iterations ("dcaplus-converging"),
// - a DCA cluster self-energy for which no lattice self-energy can be found
//   ("dca-not-terminating").

#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_sp.hpp"

#include <limits>
#include <string>

#include "dca/config/cmake_options.hpp"
#include "dca/config/threading.hpp"

#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

using namespace dca;

class LatticeMappingSpTest : public ::testing::TestWithParam<std::string> {
protected:
  using FreqDmn = func::dmn_0<phys::domains::frequency_domain>;
  using BandDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SpinDmn = func::dmn_0<phys::domains::electron_spin_domain>;
  using BandSpinDmn = func::dmn_variadic<BandDmn, SpinDmn>;

  using PointGroup = phys::domains::D4;
  using Lattice = phys::models::square_lattice<PointGroup>;
  using Model = phys::models::TightBindingModel<Lattice>;

  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, Threading, profiling::NullProfiler, Model,
                               void /*RandomNumberGenerator*/, phys::solver::CT_AUX>;
  using KClusterDmn = func::dmn_0<
      phys::domains::cluster_domain<double, Lattice::DIMENSION, phys::domains::CLUSTER,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;
  using KHostDmn = func::dmn_0<
      phys::domains::cluster_domain<double, Lattice::DIMENSION, phys::domains::LATTICE_SP,
                                    phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;

  LatticeMappingSpTest()
      : reader_(),
        sigma_cluster_("Self_Energy"),
        sigma_lattice_("Self-energy-lattice"),
        sigma_lattice_baseline_("Self-energy-lattice"),
        sigma_lattice_interpolated_("Sigma_lattice_interpolated"),
        sigma_lattice_interpolated_baseline_("Sigma_lattice_interpolated"),
        sigma_lattice_coarsegrained_("Sigma_lattice_coarsegrained"),
        sigma_lattice_coarsegrained_baseline_("Sigma_lattice_coarsegrained") {}

  static ConcurrencyType concurrency_;
  static ParametersType parameters_;

  io::HDF5Reader reader_;

  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KClusterDmn, FreqDmn>>
      sigma_cluster_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_baseline_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_interpolated_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_interpolated_baseline_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_coarsegrained_;
  func::function<std::complex<double>, func::dmn_variadic<BandSpinDmn, BandSpinDmn, KHostDmn, FreqDmn>>
      sigma_lattice_coarsegrained_baseline_;

public:
  static void SetUpTestCase() {
    parameters_.read_input_and_broadcast<io::JSONReader>(
        DCA_SOURCE_DIR "/test/integration/phys/lattice_mapping/lattice_mapping_sp_test_input.json");
    parameters_.update_model();
    parameters_.update_domains();
  }
};

LatticeMappingSpTest::ConcurrencyType LatticeMappingSpTest::concurrency_(0, nullptr);
LatticeMappingSpTest::ParametersType LatticeMappingSpTest::parameters_(
    "", LatticeMappingSpTest::concurrency_);

TEST_P(LatticeMappingSpTest, Execute) {
  reader_.open_file(DCA_SOURCE_DIR
                    "/test/integration/phys/lattice_mapping/lattice_mapping_sp_test_baseline.hdf5");
  reader_.open_group(GetParam());
  reader_.execute(sigma_cluster_);
  reader_.execute(sigma_lattice_baseline_);
  reader_.execute(sigma_lattice_interpolated_baseline_);
  reader_.execute(sigma_lattice_coarsegrained_baseline_);
  reader_.close_group();
  reader_.close_file();

  phys::latticemapping::lattice_mapping_sp<ParametersType, KClusterDmn, KHostDmn> lattice_mapping_obj(
      parameters_);

  lattice_mapping_obj.execute(sigma_cluster_, sigma_lattice_interpolated_,
                              sigma_lattice_coarsegrained_, sigma_lattice_);

  // Compate with baseline.
  for (int i = 0; i < sigma_lattice_.size(); ++i) {
    EXPECT_NEAR(sigma_lattice_baseline_(i).real(), sigma_lattice_(i).real(),
                500 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(sigma_lattice_baseline_(i).imag(), sigma_lattice_(i).imag(),
                500 * std::numeric_limits<double>::epsilon());
  }

  for (int i = 0; i < sigma_lattice_interpolated_.size(); ++i) {
    EXPECT_NEAR(sigma_lattice_interpolated_baseline_(i).real(),
                sigma_lattice_interpolated_(i).real(), 500 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(sigma_lattice_interpolated_baseline_(i).imag(),
                sigma_lattice_interpolated_(i).imag(), 500 * std::numeric_limits<double>::epsilon());
  }

  for (int i = 0; i < sigma_lattice_interpolated_.size(); ++i) {
    EXPECT_NEAR(sigma_lattice_coarsegrained_baseline_(i).real(),
                sigma_lattice_coarsegrained_(i).real(), 500 * std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(sigma_lattice_coarsegrained_baseline_(i).imag(),
                sigma_lattice_coarsegrained_(i).imag(), 500 * std::numeric_limits<double>::epsilon());
  }
}

// The last comma in the outer parentheses suppresses the gnu-zero-variadic-macro-arguments warning.
INSTANTIATE_TEST_CASE_P(InputType, LatticeMappingSpTest,
                        ::testing::Values("dcaplus-converging", "dca-not-terminating"), );
