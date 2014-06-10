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

#include "type_definitions.h"

  static const MPI_library_type           MPI_LIBRARY_TYPE  = MPI_LIBRARY;
  static const MC_integration_method_type MC_ALGORITHM_TYPE = CT_AUX;

  typedef concurrency<MPI_LIBRARY_TYPE>    concurrency_type;
  
  typedef Parameters<concurrency_type, model, MC_ALGORITHM_TYPE> parameters_type;
  
  typedef MultiOrbitalMultiSiteStructure<parameters_type, DCA_cluster_type> MOMS_type;

  //====================================================================== Create the algorithms and parameters object from the input file

  concurrency_type& concurrency(concurrency_type::get(argc,argv,"logFile"));
  concurrency << "DCA main: starting (MPI-world set up).\n\n";

  std::string stamp = get_version();
  if(concurrency.id() == concurrency.first()){
    print_version();
    provenance<interaction_type, model, MC_ALGORITHM_TYPE, MPI_LIBRARY_TYPE>::print_on_shell(QMC_INTEGRATOR_BIT);
  }

  parameters_type parameters(get_version(), concurrency);
  
  parameters.read_input_and_broadcast(file_name);

  parameters.update_model();
  parameters.update_domains();

  //====================================================================== read in functions

  MOMS_type MOMS(parameters);
  MOMS.from_JSON();

  compute_band_structure::execute(parameters, MOMS.H_LDA, MOMS.band_structure);

  
  MOMS.band_structure.print_fingerprint();

  {
    function<int   , k_DCA> I_K("I_K");
    function<double, r_DCA> distance("R");

    function<std::complex<double>, k_DCA> S_K("S_K");
    function<std::complex<double>, r_DCA> S_R("S_R");

    function<std::complex<double>, k_DCA> alpha_K("alpha_K");
    function<std::complex<double>, r_DCA> alpha_R("alpha_R");

    function<std::complex<double>, k_domain_cut_dmn_type> S_k("S_k");
    function<std::complex<double>, k_domain_cut_dmn_type> alpha_k("alpha_k");

    function<             double , r_LDA> distance_r("r");
    function<std::complex<double>, r_LDA> S_r("S_r");

    {
      I_K = -1;
      for(int k=0; k<k_DCA::dmn_size(); k++)
	for(int l=0; l<k_domain_cut_dmn_type::dmn_size(); l++)
	  if(VECTOR_OPERATIONS::L2_NORM(k_DCA::get_elements()[k], k_domain_cut_dmn_type::get_elements()[l])<1.e-4)
	    I_K(k) = l;
    }

    {
      std::vector<double> ORIGIN(2, 0.);

      distance = -1;
      for(int r=0; r<r_DCA::dmn_size(); r++)
	distance(r) = std::sqrt(DCA_r_cluster_type::find_minimal_distance(ORIGIN, r_DCA::get_elements()[r]));
    }

    for(int k=0; k<k_DCA::dmn_size(); k++)
      S_K(k) = MOMS.Sigma(0,0,k,w::dmn_size()/2);

    FT<k_DCA, r_DCA>::execute(S_K, S_R);

    {
      function<std::complex<double>, dmn_3<nu, nu, k_DCA> > Sigma;
      function<std::complex<double>, dmn_3<nu, nu, k_domain_cut_dmn_type> > Sigma_interp;

      function<std::complex<double>, dmn_3<nu, nu, k_DCA> > alpha;
      function<std::complex<double>, dmn_3<nu, nu, r_DCA> > alpha_r;
      function<std::complex<double>, dmn_3<nu, nu, k_domain_cut_dmn_type> > alpha_interp;

      double shift = 1.;

      int w_index = w::dmn_size()/2;

      memcpy(&Sigma(0), &MOMS.Sigma(0,0,0,w_index), sizeof(std::complex<double>)*k_DCA::dmn_size()*square(s::dmn_size()*b::dmn_size()));

      {	
	DCA::transform_to<DCA::ALPHA>::execute(shift, Sigma, alpha);

	wannier_interpolation<k_DCA, k_domain_cut_dmn_type>::execute(alpha, alpha_interp);

	DCA::transform_to<DCA::ALPHA>::execute_inverse(shift, Sigma_interp, alpha_interp);
      }

      FT<k_DCA, r_DCA>::execute(alpha, alpha_r);

      for(int k=0; k<k_DCA::dmn_size(); k++){
	alpha_K(k) = alpha  (0,0,k);
	alpha_R(k) = alpha_r(0,0,k);
      }

      for(int l=0; l<k_domain_cut_dmn_type::dmn_size(); l++){
	S_k(l)     = Sigma_interp(0,0,l);
	alpha_k(l) = alpha_interp(0,0,l);
      }
    }
    
    {
      std::vector<double> ORIGIN(2, 0.);
      
      distance_r = -1;
      for(int r=0; r<r_LDA::dmn_size(); r++)
	distance_r(r) = std::sqrt(LDA_r_cluster_type::find_minimal_distance(ORIGIN, r_LDA::get_elements()[r]));
    }

    {
      function<std::complex<double>, dmn_3<nu, nu, k_DCA> > Sigma;
      function<std::complex<double>, dmn_3<nu, nu, k_LDA> > Sigma_interp;

      function<std::complex<double>, dmn_3<nu, nu, k_DCA> > alpha;
      function<std::complex<double>, dmn_3<nu, nu, r_DCA> > alpha_r;
      function<std::complex<double>, dmn_3<nu, nu, k_LDA> > alpha_interp;

      function<std::complex<double>, dmn_3<nu, nu, r_LDA> > tmp;
      
      double shift = 1.;

      int w_index = w::dmn_size()/2;

      memcpy(&Sigma(0), &MOMS.Sigma(0,0,0,w_index), sizeof(std::complex<double>)*k_DCA::dmn_size()*square(s::dmn_size()*b::dmn_size()));

      {	
	DCA::transform_to<DCA::ALPHA>::execute(shift, Sigma, alpha);

	wannier_interpolation<k_DCA, k_LDA>::execute(alpha, alpha_interp);

	DCA::transform_to<DCA::ALPHA>::execute_inverse(shift, Sigma_interp, alpha_interp);
      }

      FT<k_LDA, r_LDA>::execute(Sigma_interp, tmp);
      
      for(int r=0; r<r_LDA::dmn_size(); r++)
	S_r(r) = tmp(0,0,r);
    }

    print_data::to_JSON(std::string(parameters.get_directory()+"data.json"), parameters, I_K, distance, S_K, S_R, alpha_K, alpha_R, S_k, alpha_k, distance_r, S_r);
  }

  parameters_type::profiler_type::stop(concurrency, parameters.get_profiling_file_name());
      
  concurrency << "\n\tanalysis: ending.\n";

  return 0;
}

