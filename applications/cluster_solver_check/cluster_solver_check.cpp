// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// First solves the cluster problem by exact diagonalization and then
// using QMC.
// Usage: ./cluster_solver_check input_file.json

#include <string>
#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/config/cluster_solver_check.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "phys_library/DCA+_data/moms_w_real.hpp"
#include "phys_library/DCA+_loop/DCA_loop_data.hpp"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  std::string input_file(argv[1]);

  using w = dca::func::dmn_0<frequency_domain>;
  using b = dca::func::dmn_0<electron_band_domain>;
  using s = dca::func::dmn_0<electron_spin_domain>;
  using nu = dca::func::dmn_variadic<b, s>;  // orbital-spin index
  using k_DCA =
      dca::func::dmn_0<cluster_domain<double, Lattice::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  Concurrency concurrency(argc, argv);

  Profiler::start();

  // Print some info.
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nCluster-solver-check starting.\n"
              << "MPI-world set up: " << concurrency.number_of_processors() << " processes.\n"
              << std::endl;

#ifdef DCA_WITH_CUDA
    dca::linalg::util::printInfoDevices();
#endif  // DCA_WITH_CUDA

    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

#ifdef DCA_WITH_CUDA
  dca::linalg::util::initializeMagma();
#endif  // DCA_WITH_CUDA

  // Create the parameters object from the input file.
  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  parameters.update_model();
  parameters.update_domains();

  DCA::DCA_loop_data<ParametersType> dca_loop_data;

  // Create and initialize the DCA_data objects.
  DcaData dca_data_imag(parameters);
  dca_data_imag.initialize();
  MOMS_w_real<ParametersType> dca_data_real(parameters);

  std::string data_file_ed = parameters.get_directory() + parameters.get_ED_output_file_name();
  std::string data_file_qmc = parameters.get_directory() + parameters.get_QMC_output_file_name();

  // ED solver
  EdSolver ed_solver(parameters, dca_data_imag, dca_data_real);
  ed_solver.initialize(0);
  ed_solver.execute();
  ed_solver.finalize(dca_loop_data);

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_ed(
      dca_data_imag.Sigma);

  if (concurrency.id() == concurrency.first()) {
    ed_solver.write(data_file_ed);
  }

  // QMC solver
  // The QMC solver uses the free Greens function G0 computed by the ED solver.
  // It is passed via the dca_data_imag object.
  ClusterSolver qmc_solver(parameters, dca_data_imag);
  qmc_solver.initialize(1);  // 1 = dummy iteration number
  qmc_solver.integrate();
  qmc_solver.finalize(dca_loop_data);

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_qmc(
      dca_data_imag.Sigma);

  if (concurrency.id() == concurrency.first()) {
    dca_data_imag.write(data_file_qmc);
  }

  // Check errors
  if (concurrency.id() == concurrency.first()) {
    double max_absolute_diff = 0;
    double max_relative_diff = 0;

    std::cout << "\n\nErrors\n"
              << "b1\tb2\ts1\ts2\tK\t|Sigma_ED - Sigma_QMC|_inf\t|Sigma_ED - Sigma_QMC|_inf / "
                 "|Sigma_ED|_inf"
              << std::endl;
    for (int b1 = 0; b1 < b::dmn_size(); ++b1) {
      for (int b2 = 0; b2 < b::dmn_size(); ++b2) {
        for (int s1 = 0; s1 < s::dmn_size(); ++s1) {
          for (int s2 = 0; s2 < s::dmn_size(); ++s2) {
            for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
              double ed_min_qmc_inf = 0.;
              double ed_inf = 0.;
              for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
                ed_min_qmc_inf =
                    std::max(ed_min_qmc_inf, std::abs(Sigma_ed(b1, s1, b2, s2, k_ind, w_ind) -
                                                      Sigma_qmc(b1, s1, b2, s2, k_ind, w_ind)));
                ed_inf = std::max(ed_inf, std::abs(Sigma_ed(b1, s1, b2, s2, k_ind, w_ind)));
                max_absolute_diff = std::max(max_absolute_diff, ed_min_qmc_inf);
                max_relative_diff = std::max(max_relative_diff, ed_min_qmc_inf / ed_inf);
              }
              std::cout << b1 << "\t" << b2 << "\t" << s1 << "\t" << s2 << "\t" << k_ind << "\t"
                        << ed_min_qmc_inf << "\t\t\t" << ed_min_qmc_inf / ed_inf << std::endl;
            }
          }
        }
      }
    }

    std::cout << "\n|(Sigma_ED - Sigma_QMC)|_inf          = " << max_absolute_diff
              << "\n|(Sigma_ED - Sigma_QMC)/Sigma_ED|_inf = " << max_relative_diff << std::endl;
  }

  Profiler::stop(concurrency, parameters.get_profiling_file_name());

  if (concurrency.id() == concurrency.first())
    std::cout << "\nCluster-solver-check ending.\n" << std::endl;

  return 0;
}
