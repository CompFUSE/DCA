// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// First solves the cluster problem by exact diagonalization and then
// using QMC.
// Usage: ./cluster_solver_check input_file.json

#include <string>
#include <iostream>

#include "dca/config/cluster_solver_check.hpp"
#include "dca/config/cmake_options.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/phys/dca_data/dca_data_real_freq.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  Concurrency concurrency(argc, argv);

  try {
    std::string input_file(argv[1]);

    const bool perform_statistical_test = concurrency.number_of_processors() >= 8;

    Profiler::start();

    // Print some info.
    if (concurrency.id() == concurrency.first()) {
      dca::util::GitVersion::print();
      dca::util::Modules::print();
      dca::config::CMakeOptions::print();

#ifdef DCA_WITH_CUDA
      dca::linalg::util::printInfoDevices();
#endif  // DCA_WITH_CUDA

      std::cout
          << "\n"
          << "********************************************************************************\n"
          << "**********                    Cluster Solver Check                    **********\n"
          << "********************************************************************************\n"
          << "\n"
          << "Start time : " << dca::util::print_time() << "\n"
          << "\n"
          << "MPI-world set up: " << concurrency.number_of_processors() << " processes."
          << "\n"
          << std::endl;
    }

#ifdef DCA_WITH_CUDA
    dca::linalg::util::initializeMagma();
#endif  // DCA_WITH_CUDA

    // Create the parameters object from the input file.
    ParametersType parameters(dca::util::GitVersion::string(), concurrency);
    parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
    parameters.set_dca_iterations(2);
    parameters.update_model();
    parameters.update_domains();

    dca::phys::DcaLoopData<ParametersType> dca_loop_data;

    // Create and initialize the DCA_data objects.
    DcaDataType dca_data_imag(parameters);
    dca_data_imag.initialize();
    dca::phys::DcaDataRealFreq<ParametersType> dca_data_real(parameters);

    std::string data_file_ed = parameters.get_directory() + parameters.get_filename_ed();
    std::string data_file_qmc = parameters.get_directory() + parameters.get_filename_qmc();

    // ED solver
    EdSolver ed_solver(parameters, dca_data_imag, dca_data_real);
    ed_solver.initialize(0);
    ed_solver.execute();
    ed_solver.finalize(dca_loop_data);

    const auto Sigma_ed(dca_data_imag.Sigma);
    const int tested_frequencies = 10;
    const auto G_ed(dca::math::util::cutFrequency(dca_data_imag.G_k_w, tested_frequencies));

    if (concurrency.id() == concurrency.first()) {
      ed_solver.write(data_file_ed);
    }

    // QMC solver
    // The QMC solver uses the free Greens function G0 computed by the ED solver.
    // It is passed via the dca_data_imag object.
    ClusterSolver qmc_solver(parameters, dca_data_imag);
    qmc_solver.initialize(1);  // 1 = dummy iteration number
    qmc_solver.integrate();

    // If enabled, perform statistical test.
    double p_val = -1.;
    if (perform_statistical_test) {
      auto G_qmc(dca::math::util::cutFrequency(qmc_solver.local_G_k_w(), tested_frequencies));

      using KDmn = dca::func::dmn_0<dca::phys::domains::cluster_domain<
          double, Lattice::DIMENSION, dca::phys::domains::CLUSTER,
          dca::phys::domains::MOMENTUM_SPACE, dca::phys::domains::BRILLOUIN_ZONE>>;

      dca::func::function<double, dca::math::util::CovarianceDomain<KDmn>> covariance;
      concurrency.computeCovarianceAndAvg(G_qmc, covariance);

      if (concurrency.id() == concurrency.first()) {
        dca::math::StatisticalTesting test(G_qmc, G_ed, covariance, false);

        try {
          p_val = test.computePValue(false, concurrency.number_of_processors());
          test.printInfo("statistical_test_info.txt", true);
        }
        catch (std::logic_error& err) {
          std::cerr << "Warning: " << err.what() << std::endl;
          if (test.get_dof() >= concurrency.number_of_processors())
            std::cerr << "Not enough ranks." << std::endl;
          std::cerr << "Aborting statistical test." << std::endl;
        }
      }
    }

    qmc_solver.finalize(dca_loop_data);

    if (concurrency.id() == concurrency.first()) {
      dca_data_imag.write(data_file_qmc);
    }

    // Print errors.
    if (concurrency.id() == concurrency.first()) {
      const auto& Sigma_qmc(dca_data_imag.Sigma);
      auto errors = dca::func::util::difference(Sigma_ed, Sigma_qmc);

      std::cout << "\n|(Sigma_ED - Sigma_QMC)|_1 = " << errors.l1
                << "\n|(Sigma_ED - Sigma_QMC)|_2 = " << errors.l2
                << "\n|(Sigma_ED - Sigma_QMC)|_inf = " << errors.l_inf << std::endl;

      if (p_val != -1.)
        std::cout << "\n***\nThe p-value is " << p_val << ".\n***" << std::endl;
      else if (perform_statistical_test)
        std::cout << "\n***\nStatistical test aborted.\n***" << std::endl;
    }

    Profiler::stop(concurrency, parameters.get_filename_profiling());
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
    }
  }
  catch (const std::exception& err) {
    std::cout << "Unhandled exception in main function:\n\t" << err.what();
    concurrency.abort();
  }

  return 0;
}
