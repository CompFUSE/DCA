//====================================================================
// Copyright 2013-2015 ETH Zurich.
//
// Checks the cluster solvers.
// Usage: ./check_solver inputFileName
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
  
  // Configure the calculation by selecting type definitions
  static const LIN_ALG::device_type LIN_ALG_DEVICE = LIN_ALG::CPU;

  static const CLUSTER_SOLVER_NAMES ED_CLUSTER_SOLVER_NAME = ADVANCED_ED_CLUSTER_SOLVER;
  static const CLUSTER_SOLVER_NAMES QMC_CLUSTER_SOLVER_NAME = CT_AUX_CLUSTER_SOLVER;

  static const COMP_LIB::PARALLELIZATION_LIBRARY_NAMES
    PARALLELIZATION_LIBRARY_NAME = COMP_LIB::MPI_LIBRARY;
  using concurrency_type = COMP_LIB::parallelization<PARALLELIZATION_LIBRARY_NAME>;

  using parameters_type = Parameters<concurrency_type, model, QMC_CLUSTER_SOLVER_NAME>;
  using MOMS_type = DCA_data<parameters_type>;

  using ED_solver_type = cluster_solver<ED_CLUSTER_SOLVER_NAME, LIN_ALG_DEVICE, parameters_type, MOMS_type>;
  // using QMC_solver_type = cluster_solver<QMC_CLUSTER_SOLVER_NAME, LIN_ALG_DEVICE, parameters_type, MOMS_type>;
  using posix_QMC_solver_type = posix_qmci_integrator< cluster_solver<QMC_CLUSTER_SOLVER_NAME, LIN_ALG_DEVICE, parameters_type, MOMS_type> >;

  // Create the algorithms and parameters object from the input file
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
  DCA_calculation_data  DCA_info_struct;

  MOMS_type MOMS_imag(parameters);
  MOMS_imag.initialize();

  MOMS_w_real<parameters_type> MOMS_real(parameters);

  std::string data_file_ED  = parameters.get_directory()+parameters.get_ED_output_file_name();
  std::string data_file_QMC = parameters.get_directory()+parameters.get_QMC_output_file_name();

  // ED solver
  ED_solver_type ED_obj(parameters, MOMS_imag, MOMS_real);
  ED_obj.initialize(0);
  ED_obj.execute();
  ED_obj.finalize(DCA_info_struct);

  if(concurrency.id() == concurrency.first()) {
    ED_obj.write(data_file_ED);
  }

  // QMC solver
  posix_QMC_solver_type QMC_obj(parameters, MOMS_imag);
  QMC_obj.initialize(1);
  QMC_obj.integrate();
  QMC_obj.finalize(DCA_info_struct);

  if(concurrency.id() == concurrency.first()) {
    MOMS_imag.write(data_file_QMC);
  }

  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());

  if (concurrency.id() == concurrency.last())
    std::cout << "\nCheck ending.\n" << std::endl;

  return 0;
}
