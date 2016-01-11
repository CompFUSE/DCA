//====================================================================
// Copyright 2013-2015 ETH Zurich.
//
// Description.
// Usage.
//
// Authors: Peter Staar (taa@zurich.ibm.com), IBM Research - Zurich
//          Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich
//====================================================================

#include <string>
#include <iostream>
#include <complex>
#include <cmath>

#include "version.hpp"
#include "include_files.h"
#include "type_definitions.h"

int main(int argc,char *argv[])
{
  if (argc < 2)
    {
      std::cout << "Usage: " << argv[0] << " inputFileName" << std::endl;
      return -1;
    }

  std::string file_name(argv[1]);

  //============================================================ Configure the calculation by selecting type definitions.

  static const LIN_ALG::device_type LIN_ALG_DEVICE = LIN_ALG::CPU;

  // static const DCA::CLUSTER_SOLVER_NAMES ED_CLUSTER_SOLVER_NAME  = DCA::ED_CLUSTER_SOLVER;
  static const DCA::CLUSTER_SOLVER_NAMES ED_CLUSTER_SOLVER_NAME  = DCA::ADVANCED_ED_CLUSTER_SOLVER;
  static const DCA::CLUSTER_SOLVER_NAMES QMC_CLUSTER_SOLVER_NAME = DCA::CT_AUX_CLUSTER_SOLVER;
  // static const DCA::CLUSTER_SOLVER_NAMES QMC_CLUSTER_SOLVER_NAME = DCA::SS_CT_HYB;

  static const COMP_LIB::PARALLELIZATION_LIBRARY_NAMES PARALLELIZATION_LIBRARY_NAME = COMP_LIB::MPI_LIBRARY;

  typedef COMP_LIB::parallelization<PARALLELIZATION_LIBRARY_NAME>      concurrency_type;
  typedef Parameters<concurrency_type, model, QMC_CLUSTER_SOLVER_NAME> parameters_type;

  typedef DCA::DCA_data<parameters_type> MOMS_type;

  typedef                            DCA::cluster_solver< ED_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type>          ED_solver_type;
  // typedef                            DCA::cluster_solver<QMC_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type>         QMC_solver_type;
  typedef DCA::posix_qmci_integrator<DCA::cluster_solver<QMC_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type> > posix_QMC_solver_type;
  // typedef                            DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG_DEVICE, parameters_type, MOMS_type>         HTS_solver_type;

  //====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type concurrency(argc, argv);

  if (concurrency.id() == concurrency.first())
    {
      std::cout << "DCA main: starting.\n"
                << "MPI-world set up: " << concurrency.number_of_processors() << " processes.\n" << std::endl;

      Version::print();
    }

  parameters_type::profiler_type::start();

  parameters_type parameters(Version::string(), concurrency);
  parameters.read_input_and_broadcast(file_name);
  parameters.update_model();
  parameters.update_domains();

  //====================================================================== Build the initial self energies

  DCA::DCA_calculation_data  DCA_info_struct;

  MOMS_type MOMS_imag(parameters);
  MOMS_imag.initialize();

  MOMS_w_real<parameters_type> MOMS_real(parameters);

  std::string data_file_ED  = parameters.get_directory()+parameters.get_ED_output_file_name();
  std::string data_file_QMC = parameters.get_directory()+parameters.get_QMC_output_file_name();


  ED_solver_type ED_obj(parameters, MOMS_imag, MOMS_real);

  ED_obj.initialize(0);
  ED_obj.execute();
  ED_obj.finalize(DCA_info_struct);

  if (concurrency.id() == concurrency.first())
    ED_obj.write(data_file_ED);

  FUNC_LIB::function<std::complex<double>, w> Sigma_ED;
  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    Sigma_ED(w_ind) = MOMS_imag.Sigma(0, 0, 0, w_ind);
  

  posix_QMC_solver_type QMC_obj(parameters, MOMS_imag);

  QMC_obj.initialize(1);
  QMC_obj.integrate();
  QMC_obj.finalize(DCA_info_struct);

  if(concurrency.id() == concurrency.first())
    MOMS_imag.write(data_file_QMC);

  FUNC_LIB::function<std::complex<double>, w> Sigma_QMC;
  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    Sigma_QMC(w_ind) = MOMS_imag.Sigma(0, 0, 0, w_ind);

  
  FUNC_LIB::function<double, w> Sigma_diff;
  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind)
    Sigma_diff(w_ind) = std::abs(Sigma_ED(w_ind) - Sigma_QMC(w_ind));

  const int N = 8;
  double max_abs = 0.;
  for (int w_ind = w::dmn_size()/2-N; w_ind < w::dmn_size()/2+N; ++w_ind)
    max_abs = max_abs < Sigma_diff(w_ind) ? Sigma_diff(w_ind) : max_abs;

  std::cout << "\n|Sigma_ED - Sigma_QMC|_inf = " << max_abs << std::endl;

  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());

  if (concurrency.id() == concurrency.last())
    std::cout << "\n\nDCA main: ending.\n" << std::endl;

  return 0;
}
