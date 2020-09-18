// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// System-level test for analysis.
//
// Configuration details:
// - BSE solver in DCA mode,
// - channel: particle-particle-up-down,
// - momentum transfer \vec{q}=0, frequency transfer \nu=0,
// - diagonalization of symmetric \sqrt{\chi_0}*\Gamma*\sqrt{\chi_0} matrix,
// - tight-binding model on 2D square lattice.

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

TEST(AnalysisDCAParticleParticleUpDownSymmetricTest, LeadingEigenvalues) {
  using DcaPointGroupType = phys::domains::D4;
  using LatticeType = phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = phys::models::TightBindingModel<LatticeType>;
  using ConcurrencyType = parallel::NoConcurrency;
  using ParametersType =
      phys::params::Parameters<ConcurrencyType, parallel::NoThreading, profiling::NullProfiler,
                               ModelType, void /*RngType*/, phys::solver::CT_AUX>;
  using DcaDataType = phys::DcaData<ParametersType>;
  using BseSolverType = phys::analysis::BseSolver<ParametersType, DcaDataType>;

  std::cout << "Analysis test starting: " << dca::util::print_time() << "\n" << std::endl;

  ConcurrencyType concurrency(0, nullptr);

  ParametersType parameters("", concurrency);
  parameters.read_input_and_broadcast<io::JSONReader>(
      DCA_SOURCE_DIR "/test/system-level/analysis/dca_pp_up_down_sym/input_tp.json");
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();
  dca_data.read(static_cast<std::string>(DCA_SOURCE_DIR "/test/system-level/analysis/dca_tp.hdf5"));

  BseSolverType bse_solver(parameters, dca_data);
  bse_solver.calculateSusceptibilities();

  std::cout << "\nChecking data.\n" << std::endl;

  auto& leading_eigenvalues = bse_solver.get_leading_eigenvalues();
  auto& leading_eigenvectors = bse_solver.get_leading_eigenvectors();
  auto& leading_symmetry_decomposition = bse_solver.get_leading_symmetry_decomposition();

  // Read eigenvalues, eigenvectors, and symmetry decomposition from check.hdf5.
  std::remove_const<std::remove_reference<decltype(leading_eigenvalues)>::type>::type
      leading_eigenvalues_check("leading-eigenvalues");
  std::remove_const<std::remove_reference<decltype(leading_eigenvectors)>::type>::type
      leading_eigenvectors_check("leading-eigenvectors");
  std::remove_const<std::remove_reference<decltype(leading_symmetry_decomposition)>::type>::type
      leading_symmetry_decomposition_check("leading-symmetry-decomposition");
  io::HDF5Reader reader;
  reader.open_file(DCA_SOURCE_DIR "/test/system-level/analysis/dca_pp_up_down_sym/check.hdf5");
  reader.open_group("analysis-functions");
  reader.execute(leading_eigenvalues_check);
  reader.execute(leading_eigenvectors_check);
  reader.execute(leading_symmetry_decomposition_check);
  reader.close_file();

  // Compare the computed eigenvalues, eigenvectors, and symmetry decomposition with the expected
  // result.
  for (int i = 0; i < leading_eigenvalues.size(); ++i) {
    EXPECT_NEAR(leading_eigenvalues_check(i).real(), leading_eigenvalues(i).real(), 1.e-12);
    EXPECT_NEAR(leading_eigenvalues_check(i).imag(), leading_eigenvalues(i).imag(), 1.e-12);
  }
  for (int i = 0; i < leading_eigenvectors.size(); ++i) {
    EXPECT_NEAR(std::abs(leading_eigenvectors_check(i).real()),
                std::abs(leading_eigenvectors(i).real()), 1.e-12);
    EXPECT_NEAR(std::abs(leading_eigenvectors_check(i).imag()),
                std::abs(leading_eigenvectors(i).imag()), 1.e-12);
  }
  for (int i = 0; i < leading_symmetry_decomposition.size(); ++i) {
    EXPECT_NEAR(std::abs(leading_symmetry_decomposition_check(i).real()),
                std::abs(leading_symmetry_decomposition(i).real()), 1.e-8);
    EXPECT_NEAR(std::abs(leading_symmetry_decomposition_check(i).imag()),
                std::abs(leading_symmetry_decomposition(i).imag()), 1.e-8);
  }

  std::cout << "\nWriting data.\n" << std::endl;
  bse_solver.write();

  std::cout << "\nAnalysis test ending: " << dca::util::print_time() << std::endl;
}
