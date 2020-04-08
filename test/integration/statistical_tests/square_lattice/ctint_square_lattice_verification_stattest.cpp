// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Andrei Plamada (plamada@phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Verification test of CT-AUX against a reference run

#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "test/integration/statistical_tests/square_lattice/square_lattice_setup.hpp"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(CtintVerificationTest, GreensFunction) {
  using namespace dca::testing;

  const int id = dca_test_env->concurrency.id();
  const int number_of_samples = dca_test_env->concurrency.number_of_processors();

  if (id == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType<dca::testing::CT_INT> parameters(dca::util::GitVersion::string(),
                                                  dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DcaData<dca::testing::CT_INT> data(parameters);
  data.initialize();

  // Do one QMC iteration
  ThreadedSolver<dca::testing::CT_INT> qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  // Read the portion of the Greens's function to be tested.
  using dca::func::function;
  function<double, SigmaCutDomain> G_on_node = cutFrequency(qmc_solver.local_G_k_w(), n_frequencies);

  function<double, SigmaCutDomain> G_k_w_measured(G_on_node, "G_k_w");
  dca_test_env->concurrency.sum_and_average(G_k_w_measured);

  // End of concurrent section.
  if (id == dca_test_env->concurrency.first()) {
    // read the stored reference data
    function<double, CovarianceDomain> G_k_w_covariance("G_k_w_covariance");
    function<double, SigmaCutDomain> G_k_w_expected("G_k_w");
    dca::io::HDF5Reader reader;
    reader.open_file("verification_covariance_input.hdf5");
    reader.open_group("functions");
    reader.execute(G_k_w_covariance);
    reader.execute(G_k_w_expected);
    reader.close_group();
    reader.open_group("parameters");
    int reference_n_meas;
    reader.execute("measurements_per_node", reference_n_meas);
    EXPECT_EQ(reference_n_meas, parameters.get_measurements().back());
    reader.close_file();

    dca::math::StatisticalTesting test(G_k_w_measured, G_k_w_expected, G_k_w_covariance, 1);
    double p_value = test.computePValue(true, number_of_samples);
    test.printInfo("ctint_verification_testinfo.out", true);
    double p_value_default = 0.05;
    std::cout << "\n***\nThe p-value is " << p_value << "\n***\n";
    EXPECT_LT(p_value_default, p_value);
  }

  // If many MPI ranks where used store  covariance  and mean for future testing.
  if (number_of_samples > G_k_w_measured.size()) {
    function<double, CovarianceDomain> covariance_measured("G_k_w_covariance");
    dca_test_env->concurrency.computeCovariance(G_on_node, G_k_w_measured, covariance_measured);
    if (id == dca_test_env->concurrency.last()) {
      std::cout << "\nProcessor " << id << " is writing  the covariance" << std::endl;
      dca::io::HDF5Writer writer;
      writer.open_file("verification_covariance_output.hdf5");
      writer.open_group("functions");
      writer.execute(covariance_measured);
      writer.execute(G_k_w_measured);
      writer.close_group();
      // store the number of used measurements
      writer.open_group("parameters");
      writer.execute("measurements_per_node", parameters.get_measurements().back());
      writer.execute("nodes", number_of_samples);
      writer.close_group();
      writer.close_file();
    }
  }

  // write  integrator output
  // INTERNAL: do we need it?
  //  qmc_solver.finalize();
  //  if (id == 0) {
  //    std::cout << "\nProcessor " << id << " is writing data " << std::endl;
  //    dca::io::HDF5Writer writer;
  //    writer.open_file("ctint_verification_solver_output.hdf5");
  //    writer.open_group("functions");
  //    qmc_solver.write(writer);
  //    writer.close_group();
  //    writer.close_file();
  //  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      argc, argv, dca::testing::test_directory + "square_lattice_input.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
