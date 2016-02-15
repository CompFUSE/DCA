//====================================================================
// Copyright 2013-2016 ETH Zurich.
//
// First solves the cluster problem by exact diagonalization and then
// using phtreaded CT-AUX.
// CPU version.
// Usage: ./ctaux_pthreaded_cpu_test inputFileName
//
// Authors: Peter Staar (taa@zurich.ibm.com), IBM Research - Zurich
//          Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich
//====================================================================

#include <string>
#include <iostream>

#include "gitVersion.hpp"
#include "modules.hpp"
#include "include_files.h"
#include "type_definitions.h"

using namespace DCA;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " inputFileName" << std::endl;
    return -1;
  }

  std::string file_name(argv[1]);

  // Use MPI
  using concurrency_type = COMP_LIB::parallelization<COMP_LIB::MPI_LIBRARY>;

  // Parameters and data types
  using parameters_type =
      Parameters<concurrency_type, model, CT_AUX_CLUSTER_SOLVER>;
  using MOMS_type = DCA_data<parameters_type>;

  // Cluster solvers
  using ED_solver_type =
      cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, LIN_ALG::CPU, parameters_type,
                     MOMS_type>;
  using posix_QMC_solver_type =
      posix_qmci_integrator<cluster_solver<CT_AUX_CLUSTER_SOLVER, LIN_ALG::CPU,
                                           parameters_type, MOMS_type> >;

  concurrency_type concurrency(argc, argv);

  parameters_type::profiler_type::start();

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nCheck starting.\n"
              << "MPI-world set up: " << concurrency.number_of_processors()
              << " processes.\n" << std::endl;

    GitVersion::print();
    Modules::print();
  }

  parameters_type parameters(GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast(file_name);
  parameters.update_model();
  parameters.update_domains();

  // Build the initial self energies
  DCA_calculation_data DCA_info_struct;

  MOMS_type MOMS_imag(parameters);
  MOMS_imag.initialize();

  MOMS_w_real<parameters_type> MOMS_real(parameters);

  std::string data_file_ED =
      parameters.get_directory() + parameters.get_ED_output_file_name();
  std::string data_file_QMC =
      parameters.get_directory() + parameters.get_QMC_output_file_name();

  // ED solver
  ED_solver_type ED_obj(parameters, MOMS_imag, MOMS_real);
  ED_obj.initialize(0);
  ED_obj.execute();
  ED_obj.finalize(DCA_info_struct);

  if (concurrency.id() == concurrency.first()) {
    ED_obj.write(data_file_ED);
  }

  // pthreaded CT-AUX
  posix_QMC_solver_type QMC_obj(parameters, MOMS_imag);
  QMC_obj.initialize(1);
  QMC_obj.integrate();
  QMC_obj.finalize(DCA_info_struct);

  if (concurrency.id() == concurrency.first()) {
    MOMS_imag.write(data_file_QMC);
  }

  parameters_type::profiler_type::stop(concurrency,
                                       parameters.get_profiling_file_name());

  if (concurrency.id() == concurrency.last())
    std::cout << "\nCheck ending.\n" << std::endl;

  return 0;
}
