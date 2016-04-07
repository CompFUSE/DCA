//-*-C++-*-

#include "compiler_directives.h"

#include "include_files.h"

std::string get_version()
{
  string str("DEFAULT");
  return str;
}

void print_version()
{
  string str = get_version();

  cout << "\n\n\n";
  cout << "*************************************************************************************\n";
  cout << "***                                  VERSION                                      ***\n";
  cout << "*************************************************************************************\n";
  cout << "\n\n\n";

  cout << str << endl;
}

bool file_exists(std::string file_name)
{
  std::wifstream tmp(file_name.c_str());

  if(!tmp or !tmp.good() or tmp.bad())
    {
      return false;
    }
  else
    {
      return true;
    }
}

int main(int argc,char *argv[])
{
  if (argc < 2) 
    {
      std::cout << "Usage: "<<argv[0]<<" inputFileName\n";
      return -1;
    }

  std::string file_name(argv[1]);

  //============================================================ Configure the calculation by selecting type definitions.


  /*
  static const LIN_ALG ::device_type                   LIN_ALG_DEVICE               = LIN_ALG::CPU;
  static const DCA     ::CLUSTER_SOLVER_NAMES          CLUSTER_SOLVER_NAME          = DCA::CT_AUX_CLUSTER_SOLVER;
  static const COMP_LIB::PARALLELIZATION_LIBRARY_NAMES PARALLELIZATION_LIBRARY_NAME = COMP_LIB::MPI_LIBRARY;//COMP_LIB::SERIAL_LIBRARY;//

  typedef COMP_LIB::parallelization<PARALLELIZATION_LIBRARY_NAME>  concurrency_type;
  typedef Parameters<concurrency_type, model, CLUSTER_SOLVER_NAME> parameters_type;

  //typedef MultiOrbitalMultiSiteStructure<parameters_type> MOMS_type;
  typedef DCA::DCA_data<parameters_type>          MOMS_type;

  typedef                            DCA::cluster_solver<DCA::ED_CLUSTER_SOLVER      , LIN_ALG_DEVICE, parameters_type, MOMS_type>          ED_solver_type;
  typedef                            DCA::cluster_solver<DCA::CT_AUX_CLUSTER_SOLVER  , LIN_ALG_DEVICE, parameters_type, MOMS_type>         QMC_solver_type;
  typedef DCA::posix_qmci_integrator<DCA::cluster_solver<DCA::CT_AUX_CLUSTER_SOLVER  , LIN_ALG_DEVICE, parameters_type, MOMS_type> > posix_QMC_solver_type;
  typedef                            DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG_DEVICE, parameters_type, MOMS_type>         HTS_solver_type;
  */

  static const LIN_ALG ::device_type                   LIN_ALG_DEVICE     = LIN_ALG::CPU;

  //static const DCA::CLUSTER_SOLVER_NAMES ED_CLUSTER_SOLVER_NAME  = DCA::ED_CLUSTER_SOLVER;
  static const DCA::CLUSTER_SOLVER_NAMES ED_CLUSTER_SOLVER_NAME  = DCA::ADVANCED_ED_CLUSTER_SOLVER;
  static const DCA::CLUSTER_SOLVER_NAMES QMC_CLUSTER_SOLVER_NAME = DCA::CT_AUX_CLUSTER_SOLVER;
  //static const DCA::CLUSTER_SOLVER_NAMES QMC_CLUSTER_SOLVER_NAME = DCA::SS_CT_HYB;
  
  static const COMP_LIB::PARALLELIZATION_LIBRARY_NAMES PARALLELIZATION_LIBRARY_NAME = COMP_LIB::MPI_LIBRARY;

  typedef COMP_LIB::parallelization<PARALLELIZATION_LIBRARY_NAME>      concurrency_type;
  typedef Parameters<concurrency_type, model, QMC_CLUSTER_SOLVER_NAME> parameters_type;

  typedef DCA::DCA_data<parameters_type>          MOMS_type;

  typedef                            DCA::cluster_solver< ED_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type>          ED_solver_type;
  typedef                            DCA::cluster_solver<QMC_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type>         QMC_solver_type;
  typedef DCA::posix_qmci_integrator<DCA::cluster_solver<QMC_CLUSTER_SOLVER_NAME     , LIN_ALG_DEVICE, parameters_type, MOMS_type> > posix_QMC_solver_type;
  typedef                            DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, LIN_ALG_DEVICE, parameters_type, MOMS_type>         HTS_solver_type;

  //====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type concurrency(argc, argv);

  if(concurrency.id() == concurrency.first())
    cout << "DCA main: starting (MPI-world set up).\n\n";

  parameters_type::profiler_type::start();

  std::string stamp = get_version();
  if(concurrency.id() == concurrency.first()){
    print_version();
    //provenance<interaction_type, model, MC_ALGORITHM_TYPE, MPI_LIBRARY_TYPE>::print_on_shell(QMC_INTEGRATOR_BIT);
  }

  parameters_type parameters(stamp, concurrency);

  parameters.read_input_and_broadcast(file_name);

  parameters.update_model();
  parameters.update_domains();

  //====================================================================== Build the initial self energies

  //dca::DCA_information_structure  DCA_info_struct;
  DCA::DCA_calculation_data  DCA_info_struct;

  MOMS_type MOMS_imag(parameters);
  MOMS_imag.initialize();

  MOMS_w_real<parameters_type> MOMS_real(parameters);

//   std::string data_file_ED  = parameters.get_directory()+parameters.get_output_file_name();
//   std::string data_file_CPE = parameters.get_directory()+parameters.get_spectrum_file_name();

  std::string data_file_ED  = parameters.get_directory()+parameters.get_ED_output_file_name();
  std::string data_file_CPE = parameters.get_directory()+parameters.get_CPE_output_file_name();
  std::string data_file_QMC = parameters.get_directory()+parameters.get_QMC_output_file_name();

  bool data_file_ED_exists = file_exists(data_file_ED);

  if(not data_file_ED_exists)
    {
      ED_solver_type ED_obj(parameters, MOMS_imag, MOMS_real);

      ED_obj.initialize(0);

      ED_obj.execute();

      ED_obj.finalize(DCA_info_struct);

      if(concurrency.id() == concurrency.first())
        ED_obj.write(data_file_ED);
    }

  if(parameters.do_CPE())
    {
      DCA::compute_spectrum<parameters_type, DCA::piece_wise_linear_function> CS_obj(parameters);

      if(data_file_ED_exists)
	{
	  MOMS_imag.read(data_file_ED);
          MOMS_real.read(data_file_ED);
	}

      {// reset the interacting fucntions
        MOMS_real.Sigma = 0;
        MOMS_real.G_k_w = 0;

        MOMS_real.G_k_w = 0;
        MOMS_real.G_r_w = 0;

        MOMS_real.A_w    = 0;
        MOMS_real.A_nu_w = 0;
      }

      CS_obj.execute(MOMS_imag, MOMS_real);

      if(concurrency.id() == concurrency.first())
        CS_obj.write(data_file_CPE, MOMS_imag, MOMS_real);
    }

  if(true)
    {
      posix_QMC_solver_type QMC_obj(parameters, MOMS_imag);

      QMC_obj.initialize(1);

      QMC_obj.integrate();

      QMC_obj.finalize(DCA_info_struct);

      if(concurrency.id() == concurrency.first())
        MOMS_imag.write(data_file_QMC);
    }

  /*
    {
    HTS_solver_type HTS_solver(parameters, MOMS);

    HTS_solver.initialize(2);

    HTS_solver.execute();

    HTS_solver.finalize(DCA_info_struct);

    if(concurrency.id() == concurrency.first())
    HTS_solver.write("data_HTS.json");
    }
  */


  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());

  if(concurrency.id() == concurrency.last())
    cout << "\n\nDCA main: ending. \n\n";

  return 0;
}
