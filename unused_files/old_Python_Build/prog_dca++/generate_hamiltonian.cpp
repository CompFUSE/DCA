//-*-C++-*-

#include "compiler_directives.h"

#include "include_files.h"

std::string get_version()
{
  string str(
	"commit 8e89961784846ea9d71e7525dac006d9b3870e6e\n"
	"Author: Peter Staar <staarp@itp.phys.ethz.ch>\n"
	"Date:   Tue Jul 31 15:29:54 2012 +0200\n"
	"\n"
	"compiles with -pedantic flag	2012-07-31	15:29\n");
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

int main(int argc,char *argv[])
{
  if (argc < 2) {
    std::cout << "Usage: "<<argv[0]<<" inputFileName\n";
    return -1;
  }

  std::string file_name(argv[1]);

  //============================================================ Configure the calculation by selecting type definitions.


  static const MPI_library_type            MPI_LIBRARY_TYPE    = MPI_LIBRARY;
  static const  MC_integration_method_type MC_ALGORITHM_TYPE   = CT_AUX;

  typedef concurrency<MPI_LIBRARY_TYPE>    concurrency_type;

  typedef Parameters<concurrency_type, model, MC_ALGORITHM_TYPE> parameters_type;
  
  typedef MultiOrbitalMultiSiteStructure<parameters_type, DCA_cluster_type> MOMS_type;
  
  typedef QMC::MC_integrator<MC_ALGORITHM_TYPE, parameters_type, MOMS_type> Monte_Carlo_Integrator_type;
  
  typedef dca::DCA_calculation<parameters_type, MOMS_type, concurrency_type, Monte_Carlo_Integrator_type> DCA_calculation_type;  

  //====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type& concurrency(concurrency_type::get(argc,argv,"logFile"));

  if (concurrency.id() == concurrency.first())
    cout << "DCA main: starting (MPI-world set up).\n\n";

  std::string stamp = get_version();
  if(concurrency.id() == concurrency.first()){
    print_version();
    provenance<interaction_type, model, MC_ALGORITHM_TYPE, MPI_LIBRARY_TYPE>::print_on_shell(QMC_INTEGRATOR_BIT);
  }

  parameters_type parameters(stamp, concurrency);

  parameters.read_input_and_broadcast(file_name);

  cout << "\n\n update model \n\n";
  parameters.update_model();

  cout << "\n\n update domain \n\n";
  parameters.update_domains();

  function<double              , nu_k_cut> band_structure("band-structure");

  function<int                 , nu_nu>               H_symmetry    ("H_symmetry");
  function<double              , dmn_3<nu,nu,r_DCA> > H_interactions("H_interactions");

  function<std::complex<double>, dmn_3<nu,nu,k_DCA> > H_DCA("H_DCA");
  function<std::complex<double>, dmn_3<nu,nu,k_LDA> > H_LDA("H_LDA");

  function<std::complex<double>, dmn_3<nu,nu,r_DCA> > H_r_DCA("H_r");

  {
    cout << "\n\n initialize H_symmetries \n\n";
    model::initialize_H_symmetries(H_symmetry);

    cout << "\n\n initialize H_interaction \n\n";
    model::initialize_H_interaction(H_interactions, parameters);

    cout << "\n\n initialize H_LDA \n\n";
    model::initialize_H_LDA(H_LDA, parameters);
  }

  cout << "\n\n wannier interpolation \n\n";
  wannier_interpolation<k_LDA, k_DCA>::execute(H_LDA, H_DCA);
  
  cout << "\n\n band-structure \n\n";
  compute_band_structure::execute(parameters, H_DCA, band_structure);

  cout << "\n\n wannier interpolation \n\n";
  FT<k_DCA, r_DCA>::execute(H_DCA, H_r_DCA);

  print_data::print_Hamiltionian_to_JSON(parameters.get_output_file_name(), parameters, band_structure, 
					 H_symmetry, H_interactions, H_DCA, H_r_DCA);
  
  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());

  if (concurrency.id() == concurrency.last())
    cout << "\n\nDCA main: ending. \n\n";

  return 0;
}
