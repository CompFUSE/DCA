// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// System-level test for analysis.
//
// Configuration details:
// - BSE solver in DCA mode,
// - channel: particle-particle-up-down,
// - momentum transfer \vec{q}=0, frequency transfer \nu=0,
// - diagonalization of full \Gamma*\chi_0 matrix,
// - tight-binding model on 2D square lattice.

#define DCA_ANALYSIS_TEST_WITH_FULL_DIAGONALIZATION

#include <iostream>
#include <string>
#include <utility>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_solver.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"

using namespace dca;

TEST(AnalysisDCAParticleParticleUpDownFullTest, LeadingEigenvalues) {
  using DcaPointGroupType = phys::domains::D4;
  using LatticeType = phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = phys::models::TightBindingModel<LatticeType>;
  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, parallel::NoThreading, profiling::NullProfiler,
                               ModelType, void /*RngType*/, phys::solver::CT_AUX>;
  using DcaDataType = phys::DcaData<ParametersType>;

  std::cout << "Analysis test starting.\n" << std::endl;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);
  parameters.read_input_and_broadcast<io::JSONReader>(
      DCA_SOURCE_DIR "/test/system-level/analysis/dca_pp_up_down_full/input_tp.json");
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();
  dca_data.read(static_cast<std::string>(DCA_SOURCE_DIR "/test/system-level/analysis/dca_tp.hdf5"));

  phys::analysis::BseSolver<ParametersType, DcaDataType> analysis_obj(parameters, dca_data);
  analysis_obj.calculate_susceptibilities_2();

  std::cout << "\nChecking data.\n" << std::endl;

  const static int num_lambdas = 10;
  using lambda_dmn_type = func::dmn_0<func::dmn<num_lambdas, int>>;

  func::function<std::complex<double>, lambda_dmn_type>& leading_eigenvalues =
      analysis_obj.get_leading_eigenvalues();

  // Read eigenvalues from check_data file.
  func::function<std::complex<double>, lambda_dmn_type> leading_eigenvalues_check(
      "leading-eigenvalues");
  io::HDF5Reader reader;
  reader.open_file(DCA_SOURCE_DIR "/test/system-level/analysis/dca_pp_up_down_full/check.hdf5");
  reader.open_group("analysis-functions");
  reader.execute(leading_eigenvalues_check);
  reader.close_file();

  // Compare the computed eigenvalues with the expected result.
  for (int i = 0; i < lambda_dmn_type::dmn_size(); ++i) {
    EXPECT_NEAR(leading_eigenvalues_check(i).real(), leading_eigenvalues(i).real(), 1.e-14);
    EXPECT_NEAR(leading_eigenvalues_check(i).imag(), leading_eigenvalues(i).imag(), 1.e-14);
  }

  std::cout << "\nWriting data.\n" << std::endl;
  analysis_obj.write();

  std::cout << "\nAnalysis test ending." << std::endl;
}
