// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Statistical integration test for the CT-HYB solver using  a NiO lattice.
// This is a verification test with reference data taken from a full DCA computation.
// It can run with any number of MPI ranks.

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp"
#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

// needed only if solver output is written
//#include "dca/phys/dca_loop/dca_loop_data.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;
const std::string test_directory =
    DCA_SOURCE_DIR "/test/integration/statistical_tests/real_materials/";
const std::string file_coarsegraining =
    DCA_SOURCE_DIR "/test/integration/cluster_solver/ss_ct_hyb/NiO_coarse_grained.hdf5";

using dca::func::dmn_0;
using dca::func::dmn_variadic;
using RDmn =
    dmn_0<dca::phys::domains::cluster_domain<double, 3, dca::phys::domains::CLUSTER, dca::phys::domains::REAL_SPACE,
                                             dca::phys::domains::BRILLOUIN_ZONE>>;

using SigmaDomain = dca::math::util::SigmaDomain<RDmn>;
using SigmaCutDomain = dmn_variadic<dca::math::util::details::Bdmn, RDmn,
                                    dmn_0<dca::math::util::DmnWcut>, dca::math::util::RealImagDmn>;
using CovarianceDomain = dmn_variadic<SigmaCutDomain, SigmaCutDomain>;

TEST(Ni0, GS) {
  auto* concurrency = &dca_test_env->concurrency;
  const int id = concurrency->id();
  const int n_processes = concurrency->number_of_processors();

  // print git info
  if (id == 0) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  using Model = dca::phys::models::TightBindingModel<dca::phys::models::material_lattice<
      dca::phys::models::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>>;
  using Rng = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using TestParameters =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                    Threading, dca::profiling::NullProfiler, Model,
                                    Rng, dca::phys::solver::SS_CT_HYB>;
  using Data = dca::phys::DcaData<TestParameters>;
  using ImpuritySolver = dca::phys::solver::StdThreadQmciClusterSolver<
      dca::phys::solver::SsCtHybClusterSolver<dca::linalg::CPU, TestParameters, Data>>;

  TestParameters parameters(dca::util::GitVersion::string(), *concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(test_directory + "input_NiO.json");
  // override file input for file names
  parameters.set_t_ij_file_name(test_directory + "t_ij_NiO.txt");
  parameters.set_U_ij_file_name(test_directory + "U_ij_NiO_8_lit.txt");
  parameters.update_model();
  parameters.update_domains();

  // Perform the same number of measurements per rank.
  const int meas_per_process = parameters.get_measurements();
  parameters.set_measurements(meas_per_process * dca_test_env->concurrency.number_of_processors());

  Data data(parameters);
  // initialize H only. G0 is read from file afterwards.
  data.initializeH0_and_H_i();

  // Read and broadcast the rest of the initialization from full DCA results.
  if (id == 0) {
    dca::io::HDF5Reader reader;
    reader.open_file(file_coarsegraining);
    reader.open_group("functions");
    reader.execute(data.G0_k_w);
    reader.execute(data.G0_r_t);
    reader.execute(data.G_k_w);
    reader.close_group();
    reader.open_group("additional_functions");
    reader.execute(data.Sigma_cluster);
    reader.close_group();
    reader.close_file();
  }
  dca_test_env->concurrency.broadcast(data.G0_k_w);
  dca_test_env->concurrency.broadcast(data.G0_r_t);
  dca_test_env->concurrency.broadcast(data.G_k_w);
  dca_test_env->concurrency.broadcast(data.Sigma_cluster);

  ImpuritySolver solver(parameters, data);
  solver.initialize(0);
  solver.integrate();

  using dca::func::function;
  const int n_frequencies = 10;
  function<double, SigmaCutDomain> GS_measured = dca::math::util::bandDiagonal(
      dca::math::util::cutFrequency(solver.local_GS_r_w(), n_frequencies));

  auto GS_avg(GS_measured);
  GS_avg.set_name("GS_r_w");
  concurrency->sum_and_average(GS_avg);

  // end of concurrent section
  if (id == 0) {
    // read the stored reference data
    function<double, CovarianceDomain> cov("Covariance");
    function<double, SigmaCutDomain> GS_expected("GS_r_w");
    dca::io::HDF5Reader reader;
    reader.open_file(test_directory + "NiO_covariance_input.hdf5");
    reader.open_group("functions");
    reader.execute(cov);
    reader.execute(GS_expected);
    reader.close_group();
    reader.open_group("parameters");
    int reference_n_meas;
    reader.execute("measurements_per_node", reference_n_meas);
    EXPECT_EQ(reference_n_meas, meas_per_process);
    reader.close_group();
    reader.close_file();

    dca::math::StatisticalTesting test(GS_avg, GS_expected, cov);
    double p_value = test.computePValue(true, n_processes);
    test.printInfo("NiO_testinfo.out", true);

    double p_value_default = 0.01;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }

  // if many MPI ranks were used store output covariance and mean
  if (n_processes > GS_measured.size()) {
    function<double, CovarianceDomain> cov("Covariance");
    dca_test_env->concurrency.computeCovariance(GS_measured, GS_avg, cov);
    if (id == dca_test_env->concurrency.last()) {
      std::cout << "Processor " << id << " is writing the covariance\n";
      dca::io::HDF5Writer writer;
      writer.open_file("NiO_test_covariance_output.hdf5");
      writer.open_group("functions");
      writer.execute(GS_avg);
      writer.execute(cov);
      writer.close_group();
      writer.open_group("parameters");
      writer.execute("measurements_per_node", meas_per_process);
      writer.execute("nodes", n_processes);
      writer.close_group();
      writer.close_file();
    }
  }

  // Uncomment to write integrator output.
  // dca::phys::DcaLoopData<TestParameters> loop_data;
  // solver.finalize(loop_data);
  // if (id == 0) {
  //   std::cout << "\nProcessor " << id << " is writing data " << std::endl;
  //   IO::writer<IO::HDF5> writer;
  //   writer.open_file("NiO_test_output.hdf5");
  //   solver.write(writer);
  //   writer.close_file();
  // }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(argc, argv, "");
  ::testing::AddGlobalTestEnvironment(dca_test_env);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }
  result = RUN_ALL_TESTS();
  return result;
}
