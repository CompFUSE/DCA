//-*-C++-*-

#include "compiler_directives.h"

#include "include_files.h"

std::string get_version()
{
  string str(
	"commit 48743baab7ec7ed507f59fbe4d60d39e9c70cec9\n"
	"Author: Peter Staar <staarp@itp.phys.ethz.ch>\n"
	"Date:   Fri Nov 23 16:04:52 2012 +0100\n"
	"\n"
	"lets test dlaset and 2D memcpy for the GPU-walkers !	2012-11-23	16:04\n");
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

//namespace dca {

void synchronize_devices()
{}
//}

int main(int argc,char *argv[])
{
  

  if (argc < 2) {
    std::cout << "Usage: "<<argv[0]<<" inputFileName\n";
    return -1;
  }

  std::string file_name(argv[1]);

  //============================================================ Configure the calculation by selecting type definitions.


  static const LIN_ALG::device_type       DEVICE            = LIN_ALG::CPU;
  static const MPI_library_type           MPI_LIBRARY_TYPE  = MPI_LIBRARY;
  static const MC_integration_method_type MC_ALGORITHM_TYPE = CT_AUX;

  typedef concurrency<MPI_LIBRARY_TYPE>    concurrency_type;

  typedef Parameters<concurrency_type, model, MC_ALGORITHM_TYPE> parameters_type;
  
  typedef MultiOrbitalMultiSiteStructure<parameters_type, DCA_cluster_type> MOMS_type;
  
  typedef QMC::MC_integrator<ED_CLUSTER_SOLVER, DEVICE, parameters_type, MOMS_type> ED_type;
  //typedef dca::DCA_calculation<parameters_type, MOMS_type, ED_type>                            ED_calculation_type;  

  //typedef QMC::MC_integrator<CT_AUX           , DEVICE, parameters_type, MOMS_type> QMC_type;
  typedef QMC::posix_MC_integrator<QMC::MC_integrator<CT_AUX, DEVICE, parameters_type, MOMS_type> > QMC_type; 
  
  //typedef dca::DCA_calculation<parameters_type, MOMS_type, QMC_type> DCA_calculation_type;  

  //====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type& concurrency(concurrency_type::get(argc,argv,"logFile"));

  if(concurrency.id() == concurrency.first())
    cout << "DCA main: starting (MPI-world set up).\n\n";

  parameters_type::profiler_type::start();

  std::string stamp = get_version();
  if(concurrency.id() == concurrency.first()){
    print_version();
    provenance<interaction_type, model, MC_ALGORITHM_TYPE, MPI_LIBRARY_TYPE>::print_on_shell(QMC_INTEGRATOR_BIT);
  }

  parameters_type parameters(stamp, concurrency);

  parameters.read_input_and_broadcast(file_name);

  parameters.update_model();
  parameters.update_domains();

  //====================================================================== Build the initial self energies
  
  dca::DCA_information_structure  DCA_info_struct;

  MOMS_type MOMS(parameters);
  MOMS.initialize_functions();

  {
    ED_type ED_obj(parameters, MOMS);

    ED_obj.initialize(0);

    ED_obj.integrate();

    ED_obj.finalize(DCA_info_struct);

    if(concurrency.id() == concurrency.first())
      print_data::to_JSON("data_ED.json", parameters, MOMS, ED_obj);
  }

  //if(concurrency.id() == concurrency.last())
  //print_data::to_JSON("data_ED.json", parameters, MOMS);

  {
    QMC_type QMC_obj(parameters, MOMS);

    QMC_obj.initialize(1);

    QMC_obj.integrate();

    QMC_obj.finalize(DCA_info_struct);
  }

  if(concurrency.id() == concurrency.last())
    print_data::to_JSON("data_QMC.json", parameters, MOMS);

  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());
  
  if(concurrency.id() == concurrency.last())
    cout << "\n\nDCA main: ending. \n\n";

  return 0;
}
