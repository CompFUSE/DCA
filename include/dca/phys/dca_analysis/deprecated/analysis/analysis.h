//-*-C++-*-

#ifndef ANALYSIS_H
#define ANALYSIS_H

enum MC_integration_method {CT_AUX,
                            HYBRIDIZATION,
                            HYBRIDIZATION_FULL,
                            PCM,
                            ANALYSIS,
                            ANALYSIS_INTERPOLATION,
                            ANALYSIS_COMPUTE_REDUCIBLE_VERTEX,
                            HIGH_TEMPERATURE_SERIES_SOLVER,
                            ED_CLUSTER_SOLVER};
using MC_integration_method_type = MC_integration_method;

namespace dca {

  /*! 
   * \author peter staar
   */
  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  class analysis 
  {
#include "type_definitions.h"
    
    const static int N_LAMBDA = 10;
    typedef dmn_0<dmn<N_LAMBDA, int> > lambda_dmn_type;

    const static int N_HARMONIC = 3;
    typedef dmn_0<dmn<N_HARMONIC, int> > harmonics_dmn_type;

    typedef typename parameter_type::profiler_type    profiler_t;
    typedef typename parameter_type::Concurrency_Type concurrency_t;

    typedef dmn_4<b,b,k_DCA,w_VERTEX>                   eigenvector_dmn_t;
    typedef dmn_2<eigenvector_dmn_t, eigenvector_dmn_t> matrix_dmn_t;

  public:
 
    analysis(parameter_type& parameters, MOMS_type& MOMS);
    ~analysis();

    template<class stream_type>
    void to_JSON(stream_type& ss);

//     void calculate_reducible_vertex();

    void calculate_susceptibilities();
 
  private:

    void apply_symmetries();

    void apply_particle_particle_symmetry_on_G4_k_k_w_w();

    void load_G4_b_k_w__b_k_w();
    void load_G4_0_b_k_w__b_k_w();

    void compute_Gamma_b_k_w__b_k_w();
    void compute_full_chi_0_b_k_w__b_k_w();

    void calculate_eigenvalues();

    void find_harmonic_expansion(int i, int j);

    void calculate_eigenvectors();
    void find_phase_factor();

//     void compute_chi_0();
    void compute_chi();

//     void compute_V();

//     void compute_P0();


  private:    

    parameter_type& parameters;
    MOMS_type&      MOMS;
    concurrency_t&  concurrency; 

    diagrammatic_symmetries<parameter_type> diagrammatic_symmetries_obj;

  public:

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_b_k_w__b_k_w;
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_0_b_k_w__b_k_w;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> chi_b_k_w__b_k_w;
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> chi_0_b_k_w__b_k_w;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> Gamma_b_k_w__b_k_w;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> reducible_Gamma;
    
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> full_chi_0_b_k_w__b_k_w;

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> Gamma_times_full_chi_0;

//     FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> chi;

    FUNC_LIB::function<std::complex<double>, dmn_0<dmn<1, int> > > chi_q;

    //FUNC_LIB::function<std::complex<double>, b_b__b_b> chi_q_ 0;

    //make_G4_matrix                <parameter_type, MOMS_type> make_G4_obj;
    make_G4_0_matrix              <parameter_type, MOMS_type> make_G4_0_obj;

    //compute_bubble<parameter_type, k_DCA, w_VERTEX, TRAPEZIUM_INTEGRATION> make_G4_0_CG_obj;
    compute_bubble<parameter_type, k_DCA, w_VERTEX, QUADRATURE_INTEGRATION> make_G4_0_CG_obj;

    b_b_k_DCA_w_VERTEX b_b_k_DCA_w_VERTEX_domain;

    FUNC_LIB::function<std::complex<double>, dmn_4<b,b,k_DCA,w_VERTEX> >      chi_0_function;
    FUNC_LIB::function<std::complex<double>, dmn_4<b,b,k_DCA,w_VERTEX> > full_chi_0_function;

//     FUNC_LIB::function<std::complex<double>, dmn_0<dmn<N_LAMBDA, int> > >                              leading_eigenvalues;
//     FUNC_LIB::function<std::complex<double>, dmn_2<dmn_0<dmn<N_LAMBDA, int> >, dmn_0<dmn<3, int> > > > leading_symmetries;
//     FUNC_LIB::function<std::complex<double>, dmn_2<b_b_k_DCA_w_VERTEX, dmn_0<dmn<N_LAMBDA, int> > >  > leading_eigenvectors;

    FUNC_LIB::function<std::complex<double>,       lambda_dmn_type>                        leading_eigenvalues;
    FUNC_LIB::function<std::complex<double>, dmn_2<lambda_dmn_type, harmonics_dmn_type> >  leading_symmetries;
    FUNC_LIB::function<std::complex<double>, dmn_2<lambda_dmn_type, eigenvector_dmn_t> >   leading_eigenvectors;

    FUNC_LIB::function<std::complex<double>, dmn_0<dmn<N_LAMBDA, int> > >                                       leading_V;
    FUNC_LIB::function<std::complex<double>, dmn_0<dmn<N_LAMBDA, int> > >                                       leading_P0;

    eigensystem_plan<std::complex<double>, GENERAL> eigensystem_pln;
  };

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::analysis(parameter_type& parameters_in, 
								      MOMS_type& MOMS_in):
    parameters(parameters_in),
    MOMS(MOMS_in),
    concurrency(parameters.get_concurrency()),

    diagrammatic_symmetries_obj(parameters),

    G4_b_k_w__b_k_w("G4_b_k_w__b_k_w"),
    G4_0_b_k_w__b_k_w("G4_0_b_k_w__b_k_w"),

    chi_b_k_w__b_k_w("chi_b_k_w__b_k_w"),
    chi_0_b_k_w__b_k_w("chi_0_b_k_w__b_k_w"),

    Gamma_b_k_w__b_k_w("Gamma_b_k_w__b_k_w"),
    
    reducible_Gamma("reducible_vertex"),

    full_chi_0_b_k_w__b_k_w("full_chi_0_b_k_w__b_k_w"),

    Gamma_times_full_chi_0("Gamma_times_full_chi_0"),

    //chi("chi"),
    chi_q("chi_q"),

    //make_G4_obj     (parameters_in, MOMS_in),
    make_G4_0_obj   (parameters_in, MOMS_in),
    make_G4_0_CG_obj(parameters_in),

    chi_0_function     ("chi_0_function"),
    full_chi_0_function("full_chi_0_function"),

    leading_eigenvalues("leading_eigenvalues"),
    leading_symmetries("leading_symmetries"),
    leading_eigenvectors("leading_eigenvectors"),

    leading_V("leading_V"),
    leading_P0("leading_P0"),

    eigensystem_pln(square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size())
  {
    parameters.get_output_file_name() = parameters.get_output_susceptibilities_file_name();
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::~analysis()
  {}

//   template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
//   void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::initialize()
//   {
// //     for(int i=0; i<w_VERTEX::dmn_size(); i++)
// //       for(int j=0; j<w_VERTEX_EXTENDED::dmn_size(); j++)
// // 	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
// // 	  corresponding_extended_index[i] = j;

// //     for(int j=0; j<w_VERTEX_EXTENDED::dmn_size(); j++)
// //       for(int i=0; i<w_VERTEX::dmn_size(); i++)
// //       	if(std::fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
// // 	  is_compact_frequency[j] = true;
//   }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  template<class stream_type>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::to_JSON(stream_type& ss)
  {
    leading_eigenvalues.to_JSON(ss);
    ss << ",\n";

    leading_symmetries.to_JSON(ss);
    ss << ",\n";

    chi_q.to_JSON(ss);
    ss << ",\n";

    chi_0_function.to_JSON(ss);
    ss << ",\n";

    full_chi_0_function.to_JSON(ss);
    ss << ",\n";

    leading_eigenvectors.to_JSON(ss);
 }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::calculate_susceptibilities()
  {
    if(concurrency.id() == concurrency.last())
      cout << "\n\n\n\t calculate_susceptibilities \n\n\n" << endl;

    apply_symmetries();

    make_G4_0_obj.execute(*this);

    {
      for(int w1=0; w1<w_VERTEX::dmn_size(); w1++) 
	for(int k1=0; k1<k_DCA::dmn_size(); k1++)
	  for(int l2=0; l2<b::dmn_size(); l2++)
	    for(int l1=0; l1<b::dmn_size(); l1++)
	      chi_0_function(l1,l2,k1,w1) = G4_0_b_k_w__b_k_w(l1,l2,k1,w1, l1,l2,k1,w1);
    }

    make_G4_matrix<parameter_type, MOMS_type>::execute(MOMS.G4_k_k_w_w, G4_b_k_w__b_k_w);

    compute_Gamma_b_k_w__b_k_w();

    {
      if(concurrency.id() == concurrency.last())
	cout << "\ncompute chi_0\n" << endl;
     
      if(parameters.do_DCA_plus()){

	MOMS.coarsegrain_inversion_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated, MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);    

	make_G4_0_CG_obj.execute(MOMS.H_LDA, MOMS.Sigma_lattice, full_chi_0_b_k_w__b_k_w);
      }
      else
	make_G4_0_CG_obj.execute(MOMS.H_LDA, MOMS.Sigma, full_chi_0_b_k_w__b_k_w);

      for(int w1=0; w1<w_VERTEX::dmn_size(); w1++) 
	for(int k1=0; k1<k_DCA::dmn_size(); k1++)
	  for(int l2=0; l2<b::dmn_size(); l2++)
	    for(int l1=0; l1<b::dmn_size(); l1++)
	      full_chi_0_function(l1,l2,k1,w1) = full_chi_0_b_k_w__b_k_w(l1,l2,k1,w1, l1,l2,k1,w1);
    }

//     cout << "symmetrize G4_b_k_w__b_k_w" << endl;
//     symmetrize::execute(full_chi_0_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), true);

    if(concurrency.id() == concurrency.last())
      {
	calculate_eigenvalues();

	calculate_eigenvectors();
      }

    if(false)
      {
	compute_chi();
	
	cout << std::fixed << std::setprecision(6) 
	     << "\n\n\t chi_q = {" << 1./parameters.get_beta() << ", " << real(chi_q(0)) << "}\n\n\n";
      }
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::apply_symmetries()
  {
    if(concurrency.id() == concurrency.last())
      cout << __FUNCTION__ << endl;
    
    symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

    symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);

//     cout << "symmetrize MOMS.G4_k_k_w_w" << endl;
//     symmetrize::execute(MOMS.G4_k_k_w_w, MOMS.H_symmetry, parameters.get_q_vector(), false);

//     if(parameters.get_vertex_measurement_type() == PARTICLE_PARTICLE_UP_DOWN)
//       apply_particle_particle_symmetry_on_G4_k_k_w_w();
  }
    
  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_Gamma_b_k_w__b_k_w()
  {
    if(concurrency.id() == concurrency.last())
      cout << "\n\n\n\t compute_Gamma_b_k_w__b_k_w \n\n\n" << endl;

    int N = b_b_k_DCA_w_VERTEX_domain.get_size();

//     diagrammatic_symmetries_obj.execute(G4_b_k_w__b_k_w);
//     diagrammatic_symmetries_obj.execute(G4_0_b_k_w__b_k_w);

    {
//       cout << "symmetrize G4_b_k_w__b_k_w" << endl;
//       symmetrize::execute(G4_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), false);

      invert_plan<std::complex<double> > invert_pln(N);
      
      memcpy(&invert_pln.Matrix[0], &G4_b_k_w__b_k_w(0), sizeof(std::complex<double>)*N*N);
      invert_pln.execute_plan();
      memcpy(&chi_b_k_w__b_k_w(0), invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);

//       cout << "symmetrize chi_b_k_w__b_k_w" << endl;
//       symmetrize::execute(chi_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), true);
    }

    {
//       cout << "symmetrize G4_0_b_k_w__b_k_w" << endl;
//       symmetrize::execute(G4_0_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), true);

      invert_plan<std::complex<double> > invert_pln(N);
      
      memcpy(&invert_pln.Matrix[0], &G4_0_b_k_w__b_k_w(0), sizeof(std::complex<double>)*N*N);
      invert_pln.execute_plan();
      memcpy(&chi_0_b_k_w__b_k_w(0), invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);

//       cout << "symmetrize chi_0_b_k_w__b_k_w" << endl;
//       symmetrize::execute(chi_0_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), true);
    }

    for(int i=0; i<N*N; i++)
      Gamma_b_k_w__b_k_w(i) = chi_0_b_k_w__b_k_w(i)-chi_b_k_w__b_k_w(i);

//     cout << "symmetrize Gamma_b_k_w__b_k_w" << endl;
//     symmetrize::execute(Gamma_b_k_w__b_k_w, MOMS.H_symmetry, parameters.get_q_vector(), true);
//     diagrammatic_symmetries_obj.execute(Gamma_b_k_w__b_k_w);
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::calculate_eigenvalues()
  {
    cout << __FUNCTION__ << endl;

    cout << scientific;
    cout.precision(6);

    int matrix_dim = b_b_k_DCA_w_VERTEX_domain.get_size();

    {// \Chi_0 * \Gamma --> Gamma_times_full_chi_0(0);
      gemm_plan<std::complex<double> > gemm_pln(b_b_k_DCA_w_VERTEX_domain.get_size());
      
      gemm_pln.A = &Gamma_b_k_w__b_k_w(0);
      gemm_pln.B = &full_chi_0_b_k_w__b_k_w(0);
      gemm_pln.C = &Gamma_times_full_chi_0(0);
      
      gemm_pln.execute_plan();
      
//       cout << "symmetrize Gamma_times_full_chi_0" << endl;
//       symmetrize::execute(Gamma_times_full_chi_0, MOMS.H_symmetry, parameters.get_q_vector(), true);
    }

    std::vector<std::pair<std::complex<double>, int> > eigenvals;
    { 
      memcpy(&eigensystem_pln.A[0], &Gamma_times_full_chi_0(0), sizeof(std::complex<double>)*Gamma_times_full_chi_0.size()); 

      eigensystem_pln.execute_plan();
      
      for(int i=0; i<matrix_dim; i++)
	eigenvals.push_back(std::pair<std::complex<double>, int>(eigensystem_pln.W[i], i));
    }

    { // check imaginary part ..
      stable_sort(eigenvals.begin(), eigenvals.end(), &complex_less_pairs);
      
      cout << "\n\n max Imag[eigenvalues] :: \n\n";
      size_t i=0;
      while(i<eigenvals.size()-1 && real(eigenvals[i].first) > 0.95 ){
	cout << "\t\t " << i << "\t --> " << eigenvals[i].first << endl;
	i++;
      }
      cout << "\n\n";
    }
    
    { // write down leading eigenvalues ...
      cout << "\n\n sorted eigenvalues :: \n\n";
      
      stable_sort(eigenvals.begin(), eigenvals.end(), &susceptibility_less_pairs);

      for (int i=0; i<N_LAMBDA; i++)
	{
	  leading_eigenvalues(i) = eigenvals[eigenvals.size()-1-i].first;
	    
	  find_harmonic_expansion(i, eigenvals[eigenvals.size()-1-i].second);

	  cout.precision(6);
	  cout << "\t ---> (leading) j=" << i 
	       << "\t sval = "           << sqrt(square(1.0 - real(eigenvals[eigenvals.size()-1-i].first)) + square(imag(eigenvals[eigenvals.size()-1-i].first)))
	       << "\t eigenval = "       << real(eigenvals[eigenvals.size()-1-i].first) << " ; " << std::fabs(imag(eigenvals[eigenvals.size()-1-i].first));

	  cout << "\t|\t";
	  for(int psi=0; psi<3; psi++)
	    cout << "\t" << abs(leading_symmetries(i, psi));
	  cout << "\n";
	}
      cout << "\n\n";
    }
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::find_harmonic_expansion(int i, int leading_lambda_index)
  {
    int MATRIX_DIM = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();
      
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX > eigenvector;
    for(int l=0; l<MATRIX_DIM; l++)
      eigenvector(l) = eigensystem_pln.VR[l + leading_lambda_index*MATRIX_DIM];
    
    {// s-wave
      std::complex<double> coef=0;
      std::complex<double> norm_psi=0;
      std::complex<double> norm_phi=0;

      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
	  
	std::complex<double> psi = 1.;

	coef     += psi* eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2);
	norm_psi += psi*psi;
	norm_phi += (eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)*std::conj(eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)));
      }

      leading_symmetries(i, 0) = coef/sqrt(norm_phi)/sqrt(norm_psi);
    }

    {// p-wave
      std::complex<double> coef=0;
      std::complex<double> norm_psi=0;
      std::complex<double> norm_phi=0;

      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){

// 	double alpha_x = DCA_cluster_type::get_r_basis()[0][0];
// 	double alpha_y = DCA_cluster_type::get_r_basis()[1][1];
	double alpha_x = r_DCA::parameter_type::get_basis_vectors()[0][0];
	double alpha_y = r_DCA::parameter_type::get_basis_vectors()[1][1];

	std::complex<double> psi = (cos(alpha_x*k_DCA::get_elements()[k_ind][0]) + cos(alpha_y*k_DCA::get_elements()[k_ind][1]));

	coef += psi * eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2);
	norm_psi += psi*psi;
	norm_phi += (eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)*std::conj(eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)));
      }

      leading_symmetries(i, 1) = coef/sqrt(norm_phi)/sqrt(norm_psi);
    }

    {// d-wave
      std::complex<double> coef=0;
      std::complex<double> norm_psi=0;
      std::complex<double> norm_phi=0;

//       double alpha_x = DCA_cluster_type::get_r_basis()[0][0];
//       double alpha_y = DCA_cluster_type::get_r_basis()[1][1];
      double alpha_x = r_DCA::parameter_type::get_basis_vectors()[0][0];
      double alpha_y = r_DCA::parameter_type::get_basis_vectors()[1][1];

      for(int k_ind=0; k_ind<k_DCA::dmn_size(); k_ind++){
	  
	std::complex<double> psi = (cos(alpha_x*k_DCA::get_elements()[k_ind][0]) - cos(alpha_y*k_DCA::get_elements()[k_ind][1]));

	coef += psi * eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2);
	norm_psi += psi*psi;
	norm_phi += (eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)*std::conj(eigenvector(0,0,k_ind, w_VERTEX::dmn_size()/2)));
      }

      leading_symmetries(i, 2) = coef/sqrt(norm_phi)/sqrt(norm_psi);
    }
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::calculate_eigenvectors()
  {
    int MATRIX_DIM = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();
      
    for(int i=0; i<N_LAMBDA; i++)
      {
	int leading_lambda_index=0;
	for(leading_lambda_index=0; leading_lambda_index<MATRIX_DIM; leading_lambda_index++)
	  if(abs(eigensystem_pln.W[leading_lambda_index]-leading_eigenvalues(i))<1.e-10)
	    break;
	
	assert(abs(eigensystem_pln.W[leading_lambda_index]-leading_eigenvalues(i))<1.e-10);
	
	for(int l=0; l<MATRIX_DIM; l++)
	  leading_eigenvectors(i,l) = eigensystem_pln.VR[l + leading_lambda_index*MATRIX_DIM];
      }

    find_phase_factor();
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::find_phase_factor()
  {
    int MATRIX_DIM = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();

    for(int i=0; i<N_LAMBDA; i++)
      {
	double N=1.e4;
	
	std::complex<double> alpha_min=0;
	double norm = 1.e6;

	for(int l=0; l<N; l++)
	  {
	    std::complex<double> alpha = std::complex<double>(cos(2.*M_PI*l/N), sin(2.*M_PI*l/N)); 
	    
	    double result=0;
	    
	    for(int w1=0; w1<w_VERTEX::dmn_size()/2; w1++)
	      for(int K1=0; K1<k_DCA::dmn_size(); K1++)
		for(int m1=0; m1<b::dmn_size(); m1++)
		  for(int n1=0; n1<b::dmn_size(); n1++)
		    result += abs(alpha*leading_eigenvectors(i, n1,m1,K1,w1)
				  - conj(alpha*leading_eigenvectors(i, n1,m1,K1,w_VERTEX::dmn_size()-1-w1)));
	  
	    if(result < norm){
	      norm = result;
	      alpha_min = alpha;
	    }
	  }
	
// 	if(imag(alpha_min*leading_eigenvectors(i, 0,0,0,w_VERTEX::dmn_size()/2))>1.e-6)
// 	  alpha_min *= -1.;

	for(int l=0; l<MATRIX_DIM; l++)
	  leading_eigenvectors(i, l) *= alpha_min;

// 	double max = 0;
// 	std::complex<double> alpha=0;

// 	for(int K1=0; K1<k_DCA::dmn_size(); K1++)
// 	  for(int m1=0; m1<b::dmn_size(); m1++)
// 	    for(int n1=0; n1<b::dmn_size(); n1++)
// 	      if(abs(leading_eigenvectors(i, n1,m1,K1,w_VERTEX::dmn_size()/2))>max){
// 		max   = abs(leading_eigenvectors(i, n1,m1,K1,w_VERTEX::dmn_size()/2));
// 		alpha =     leading_eigenvectors(i, n1,m1,K1,w_VERTEX::dmn_size()/2);
// 	      }
	
// 	for(int l=0; l<MATRIX_DIM; l++)
// 	  leading_eigenvectors(i, l) /= alpha;
      }
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_chi()
  {
    cout << __FUNCTION__ << endl;

    // eqn 20-21 in PRB 64 195130
    // \chi(Q,T) = \frac{1}{(\beta*N_c)^2} \sum_{K1,K2} \chi[Q,K1,K2]
    //           ===> \chi[Q,K1,K2] = [1-\chi^0 * \Gamma]^{-1} * \chi^0[Q,K1,K2]

    int N = b_b_k_DCA_w_VERTEX_domain.get_size();

    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> chi("chi");
    FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> denominator("temporary");

    {// \Chi_0 * \Gamma --> Gamma_times_full_chi_0(0);
      gemm_plan<std::complex<double> > gemm_pln(N);
      
      gemm_pln.A = &full_chi_0_b_k_w__b_k_w(0);
      gemm_pln.B = &Gamma_b_k_w__b_k_w(0);
      gemm_pln.C = &denominator(0);
      
      gemm_pln.execute_plan();
    }

    cout << "\t compute denominator \n";

    eigensystem_plan<std::complex<double>, GENERAL> eigensystem_pln(N,'V','V');
    { 
      memcpy(&eigensystem_pln.A[0], &denominator(0), sizeof(std::complex<double>)*denominator.size()); 
      eigensystem_pln.execute_plan();

      {
	invert_plan<std::complex<double> > invert_pln(N);
	
	memcpy(invert_pln.Matrix , eigensystem_pln.VR        , sizeof(std::complex<double>)*N*N);
	invert_pln.execute_plan();
	memcpy(eigensystem_pln.VL, invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);
      }

      {
	for(int j=0; j<N; j++)
	  for(int i=0; i<N; i++)
	    eigensystem_pln.VL[i+j*N] = (1./(1.-eigensystem_pln.W[i]))*eigensystem_pln.VL[i+j*N];

	gemm_plan<std::complex<double> > gemm_pln(N);
	gemm_pln.A = eigensystem_pln.VR;
	gemm_pln.B = eigensystem_pln.VL;
	gemm_pln.C = &denominator(0,0);

	gemm_pln.execute_plan();
      }

      cout << "\t compute chi_k_k' \n";

      {
	gemm_plan<std::complex<double> > gemm_pln(N);
	gemm_pln.A = &denominator(0,0);
	gemm_pln.B = &full_chi_0_b_k_w__b_k_w(0,0);
	gemm_pln.C = &chi(0,0);

	gemm_pln.execute_plan();
      }  
      
      for(int i=0; i < chi.size(); i++) 
	chi_q(0) += chi(i);
    
      chi_q(0) *= 1./( parameters.get_beta() * k_DCA::dmn_size() );
    }
  }
}


#endif
