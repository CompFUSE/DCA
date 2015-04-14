//-*-C++-*-

#ifndef ANALYSIS_H
#define ANALYSIS_H

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
// // 	if(fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
// // 	  corresponding_extended_index[i] = j;

// //     for(int j=0; j<w_VERTEX_EXTENDED::dmn_size(); j++)
// //       for(int i=0; i<w_VERTEX::dmn_size(); i++)
// //       	if(fabs(w_VERTEX::parameter_type::get_elements()[i]-w_VERTEX_EXTENDED::parameter_type::get_elements()[j])<1.e-6)
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
     
      if(parameters.use_interpolated_Self_energy()){

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

//     if(parameters.get_vertex_measurement_type() == PARTICLE_PARTICLE_SUPERCONDUCTING)
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
	       << "\t eigenval = "       << real(eigenvals[eigenvals.size()-1-i].first) << " ; " << fabs(imag(eigenvals[eigenvals.size()-1-i].first));

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



/*
  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_V()
  {
    cout << __FUNCTION__ << endl;

    int MATRIX_DIM = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();
    
    for(int i=0; i<N_LAMBDA; i++){

      leading_V(i) = 0;
      std::complex<double> renorm = 0;

      for(int l1=0; l1<MATRIX_DIM; l1++){
	renorm  += std::conj(leading_eigenvectors(l1, i))*leading_eigenvectors(l1, i);

	for(int l2=0; l2<MATRIX_DIM; l2++)
	  leading_V(i) += std::conj(leading_eigenvectors(l1, i))*Gamma_b_k_w__b_k_w(l1,l2)*leading_eigenvectors(l2, i);
      }
      
      leading_V(i) /= renorm;
    }
  }

  template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
  void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_P0()
  {
    cout << __FUNCTION__ << endl;

    int MATRIX_DIM = square(b::dmn_size())*k_DCA::dmn_size()*w_VERTEX::dmn_size();

    double renorm = -1.;

    for(int i=0; i<N_LAMBDA; i++){

      leading_P0(i) = 0;

      for(int l1=0; l1<MATRIX_DIM; l1++)
	for(int l2=0; l2<MATRIX_DIM; l2++)
	  leading_P0(i) += std::conj(leading_eigenvectors(l1, i))*full_chi_0_b_k_w__b_k_w(l1,l2)*leading_eigenvectors(l2, i);

      leading_P0(i) *= renorm;
    }
  }
 */




















//   template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
//   void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_chi_0()
//   {
//     cout << __FUNCTION__ << endl;

//     for(int w2=0; w2<w_VERTEX::dmn_size(); w2++) 
//       for(int k2=0; k2<k_DCA::dmn_size(); k2++) 
// 	for(int l3=0; l3<b::dmn_size(); l3++)
// 	  for(int l4=0; l4<b::dmn_size(); l4++)

// 	    for(int w1=0; w1<w_VERTEX::dmn_size(); w1++) 
// 	      for(int k1=0; k1<k_DCA::dmn_size(); k1++)
// 		for(int l1=0; l1<b::dmn_size(); l1++)
// 		  for(int l2=0; l2<b::dmn_size(); l2++)
// 		    chi_q_0(l1,l2, l3,l4) += full_chi_0_b_k_w__b_k_w(l1,l2,k1,w1, l3,l4,k2,w2);

//     chi_q_0 /= (parameters.get_beta()*DCA_cluster_type::get_cluster_size());
//   }

//   template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
//   void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_chi()
//   {
//     cout << __FUNCTION__ << endl;

//     // eqn 20-21 in PRB 64 195130
//     // \chi(Q,T) = \frac{1}{(\beta*N_c)^2} \sum_{K1,K2} \chi[Q,K1,K2]
//     //           ===> \chi[Q,K1,K2] = [1-\chi^0 * \Gamma]^{-1} * \chi^0[Q,K1,K2]

//     int N = b_b_k_DCA_w_VERTEX_domain.get_size();

//     FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> denominator("temporary");

//     {// \Chi_0 * \Gamma --> Gamma_times_full_chi_0(0);
//       gemm_plan<std::complex<double> > gemm_pln(N);
      
//       gemm_pln.A = &full_chi_0_b_k_w__b_k_w(0);
//       gemm_pln.B = &Gamma_b_k_w__b_k_w(0);
//       gemm_pln.C = &denominator(0);
      
//       gemm_pln.execute_plan();
//     }

//     cout << "\t compute denominator \n";

//     eigensystem_plan<std::complex<double>, GENERAL> eigensystem_pln(N,'V','V');
//     { 
//       memcpy(&eigensystem_pln.A[0], &denominator(0), sizeof(std::complex<double>)*denominator.size()); 
//       eigensystem_pln.execute_plan();

//       {
// 	invert_plan<std::complex<double> > invert_pln(N);
	
// 	memcpy(invert_pln.Matrix , eigensystem_pln.VR        , sizeof(std::complex<double>)*N*N);
// 	invert_pln.execute_plan();
// 	memcpy(eigensystem_pln.VL, invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);
//       }

//       {
// 	for(int j=0; j<N; j++)
// 	  for(int i=0; i<N; i++)
// 	    eigensystem_pln.VL[i+j*N] = (1./(1.-eigensystem_pln.W[i]))*eigensystem_pln.VL[i+j*N];

// 	gemm_plan<std::complex<double> > gemm_pln(N);
// 	gemm_pln.A = eigensystem_pln.VR;
// 	gemm_pln.B = eigensystem_pln.VL;
// 	gemm_pln.C = &denominator(0,0);

// 	gemm_pln.execute_plan();
//       }

//       cout << "\t compute chi_k_k' \n";

//       {
// 	gemm_plan<std::complex<double> > gemm_pln(N);
// 	gemm_pln.A = &denominator(0,0);
// 	gemm_pln.B = &full_chi_0_b_k_w__b_k_w(0,0);
// 	gemm_pln.C = &chi(0,0);

// 	gemm_pln.execute_plan();
//       }  
      
//       for(int w2=0; w2<w_VERTEX::dmn_size(); w2++) 
// 	for(int k2=0; k2<k_DCA::dmn_size(); k2++) 
// 	  for(int l3=0; l3<b::dmn_size(); l3++)
// 	    for(int l4=0; l4<b::dmn_size(); l4++)
	      
// 	      for(int w1=0; w1<w_VERTEX::dmn_size(); w1++) 
// 		for(int k1=0; k1<k_DCA::dmn_size(); k1++)
// 		  for(int l1=0; l1<b::dmn_size(); l1++)
// 		    for(int l2=0; l2<b::dmn_size(); l2++)
// 		      chi_q(l1,l2, l3,l4) += chi(l1,l2,k1,w1, l3,l4,k2,w2);
      
//       chi_q /= (parameters.get_beta()*DCA_cluster_type::get_cluster_size());

// //       std::complex<double> result = 0.;
// //       for(int i=0; i < chi.size(); i++) 
// // 	result += chi(i);
    
// //       result *= 1/( parameters.get_beta() * DCA_cluster_type::get_cluster_size() );
    
// //       return result;
//     }
//   }








// find_phase_factor()
//   {
//     cout << __FUNCTION__ << endl;

//     // \alpha \phi(\varpi) == conj(\alpha \phi(-\varpi))

//     double N=1.e5;

//     std::complex<double> alpha_min=0;
//     double norm = 1.e6;

//     int N_w = w_VERTEX_EXTENDED::dmn_size()-1;

//     for(int l=0; l<N; l++)
//       {
// 	std::complex<double> alpha=0;
// 	real(alpha) = cos(2.*M_PI*l/N); 
// 	imag(alpha) = sin(2.*M_PI*l/N); 

// 	double result=0;

// 	for(int n1=0; n1<b::dmn_size(); n1++)
// 	  for(int m1=0; m1<b::dmn_size(); m1++)
// 	    for(int K1=0; K1<DCA_k_cluster_type::get_size(); K1++)
// 	      for(int w1=0; w1<w_VERTEX_EXTENDED::dmn_size()/2; w1++)
// 		result += abs(alpha*leading_eigenvector(n1,m1,K1,w1)-conj(alpha*leading_eigenvector(n1,m1,K1,N_w-w1)));

// 	if(result < norm){
// 	  norm = result;
// 	  alpha_min = alpha;
// 	}
//       }

//     if(imag(alpha_min*leading_eigenvector(0,0,0,w_VERTEX_EXTENDED::dmn_size()/2))>1.e-6)
//       alpha_min *= -1.;

//     for(int l=0; l<leading_eigenvector.size(); l++)
//       leading_eigenvector(l) *= alpha_min;
//   }















//   template<class parameter_type, class MOMS_type, MC_integration_method_type Monte_Carlo_solver_t>
//   void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::calculate_reducible_vertex()
//   {
//     cout << __FUNCTION__ << endl;

//     double renorm = 1./(parameters.get_beta()*DCA_cluster_type::get_cluster_size());
//     int N = b_b_k_DCA_w_VERTEX_domain.get_size();

//     apply_symmetries();

//     {
//       make_G4_0_obj.execute(*this);
//       //make_G4_0_obj.execute(MOMS.G_k_w, G4_0_b_k_w__b_k_w);
//     }

//     //make_G4_obj.execute(*this);
//     make_G4_matrix<parameter_type, MOMS_type>::execute(MOMS.G4_k_k_w_w, G4_b_k_w__b_k_w);

//     G4_b_k_w__b_k_w   *= renorm;
//     G4_0_b_k_w__b_k_w *= renorm;

//     G4_b_k_w__b_k_w -= G4_0_b_k_w__b_k_w;

//     FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> inverted_G4_0_b_k_w__b_k_w("inverted_G4_0_b_k_w__b_k_w");
//     {
//       invert_plan<std::complex<double> > invert_pln(N);
      
//       memcpy(&invert_pln.Matrix[0], &G4_0_b_k_w__b_k_w(0), sizeof(std::complex<double>)*N*N);
//       invert_pln.execute_plan();
//       memcpy(&inverted_G4_0_b_k_w__b_k_w(0), invert_pln.inverted_matrix, sizeof(std::complex<double>)*N*N);
//       symmetrize::execute(inverted_G4_0_b_k_w__b_k_w, MOMS.H_symmetry, true);
//     }

//     FUNC_LIB::function<std::complex<double>, b_b_k_DCA_w_VERTEX__b_b_k_DCA_w_VERTEX> G4_min_full_chi_0__times_inverted_full_chi_0("tmp");

//     {
//       gemm_plan<std::complex<double> > gemm_pln(N);
      
//       gemm_pln.A = &G4_b_k_w__b_k_w(0);
//       gemm_pln.B = &inverted_G4_0_b_k_w__b_k_w(0);
//       gemm_pln.C = &G4_min_full_chi_0__times_inverted_full_chi_0(0);
      
//       gemm_pln.execute_plan();
//     }

//     {
//       gemm_plan<std::complex<double> > gemm_pln(N);
      
//       gemm_pln.A = &inverted_G4_0_b_k_w__b_k_w(0);
//       gemm_pln.B = &G4_min_full_chi_0__times_inverted_full_chi_0(0);
//       gemm_pln.C = &reducible_Gamma(0);
      
//       gemm_pln.execute_plan();
//     }
//   }





















//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type, Monte_Carlo_solver_t>::compute_full_chi_0_b_k_w__b_k_w()
//   {
//     cout << __FUNCTION__ << endl;

//     int                 W     = parameters.get_w_channel();
//     //std::vector<double> Q     = DCA_k_cluster_type::get_elements()[parameters.get_q_channel()];
//     double chemical_potential = parameters.get_chemical_potential();

//     const static double               Nb_interpolation = 16;
//     std::vector<std::vector<double> > centered_mesh    = Mesh<DCA_k_cluster_type>::execute(int(Nb_interpolation));
//     std::vector<std::vector<double> > mesh             = centered_mesh;

//     double integration_factor = get_integration_factor()/double(centered_mesh.size());

//     int matrix_size = 2*2*b::dmn_size()*b::dmn_size();
//     int matrix_dim  = 2*b::dmn_size();
    
//     std::complex<double>* tmp_matrix            = new std::complex<double>[matrix_size];
//     std::complex<double>* Sigma_matrix          = new std::complex<double>[matrix_size];
//     std::complex<double>* H_LDA_matrix_k        = new std::complex<double>[matrix_size];
//     std::complex<double>* H_LDA_matrix_k_accent = new std::complex<double>[matrix_size];
//     std::complex<double>* G_k                   = new std::complex<double>[matrix_size];
//     std::complex<double>* G_k_accent            = new std::complex<double>[matrix_size];
    
//     wannier_interpolation::mesh_k::get_size() = int(mesh.size());
//     wannier_interpolation WIP_k       (MOMS.H_LDA);
//     wannier_interpolation WIP_k_accent(MOMS.H_LDA);

//     invert_plan<std::complex<double> > invert_pln(matrix_dim);

//     for(int K_ind=0; K_ind<DCA_k_cluster_type::get_size(); K_ind++)
//       {
// 	std::vector<double> K = DCA_k_cluster_type::get_elements()[K_ind];
// 	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K);
// 	wannier_interpolation::H_k_interpolated_type& H_k        = WIP_k.execute(mesh);

// 	int                 K_accent_ind = get_k_accent(K_ind).first;
// 	std::vector<double> K_accent     = get_k_accent(K_ind).second;
// 	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K_accent);
// 	wannier_interpolation::H_k_interpolated_type& H_k_accent = WIP_k_accent.execute(mesh);

// 	// \sum_{k \in K_{k}}
// 	for(int k_ind=0; k_ind<int(mesh.size()); k_ind++)
// 	  {
// 	    memcpy(H_LDA_matrix_k       , &H_k       (0,0,k_ind), sizeof(std::complex<double>)*matrix_size);
// 	    memcpy(H_LDA_matrix_k_accent, &H_k_accent(0,0,k_ind), sizeof(std::complex<double>)*matrix_size);

// 	    for(int w_vertex_index=0; w_vertex_index<w_VERTEX_EXTENDED::dmn_size(); w_vertex_index++)
// 	      {
// 		double w  = w_VERTEX_EXTENDED::parameter_type::get_elements()[w_vertex_index];
// 		int w_ind = w_vertex_index - w_VERTEX_EXTENDED::dmn_size()/2 + w::dmn_size()/2;
		
// 		assert( fabs(w -  w::parameter_type::get_elements()[w_ind]) < 1.e-6);

// 		{// G(k)
// 		  memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_ind, w_ind), sizeof(std::complex<double>)*matrix_size);

// 		  {// G^-1 =  -(H_k + Sigma) + i*w + mu
// 		    for(int index=0; index < matrix_size; index++)
// 		      tmp_matrix[index] = -(H_LDA_matrix_k[index] + Sigma_matrix[index]);  
		    
// 		    for(int nu=0; nu<matrix_dim; nu++)
// 		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w);
// 		  }
		  
// 		  {// G(k)^-1 = (-(H(k)+Sigma(k)) + i*w) ==> G(k) = (-(H(k)+Sigma(k)) + i*w)^-1 
// 		    memcpy(invert_pln.Matrix, tmp_matrix                    , sizeof(std::complex<double>)*matrix_size);
// 		    invert_pln.execute_plan();
// 		    memcpy(&G_k[0]          , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
// 		  }
// 		}

// 		{// G(k+Q)
// 		  std::pair<int, double> w_accent = get_w_accent(w_ind, W);
		  
// 		  memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_accent_ind, w_accent.first), sizeof(std::complex<double>)*matrix_size);

// 		  {// G^-1 =  -(H_k_accent + Sigma) + i*w + mu
// 		    for(int index=0; index < matrix_size; index++)
// 		      tmp_matrix[index] = -(H_LDA_matrix_k_accent[index] + Sigma_matrix[index]);  
		    
// 		    for(int nu=0; nu<matrix_dim; nu++)
// 		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w_accent.second);
// 		  }
		  
// 		  {// G(k+Q)^-1 = (-(H(k+Q)+Sigma(k+Q)) + i*w) ==> G(k+Q) = (-(H(k+Q)+Sigma(k+Q)) + i*w)^-1 
// 		    memcpy(invert_pln.Matrix, tmp_matrix                     , sizeof(std::complex<double>)*matrix_size);
// 		    invert_pln.execute_plan();
// 		    memcpy(&G_k_accent[0]    , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
// 		  }
// 		}

// 		for(int b1=0; b1<b::dmn_size(); b1++){
// 		  for(int b2=0; b2<b::dmn_size(); b2++){

// 		    std::complex<double> c(0,0);

// 		    switch(parameters.get_vertex_measurement_type()) 
// 		      {
// 		      case PARTICLE_HOLE_MAGNETIC:
// 			for(int b3=0; b3<b::dmn_size(); b3++)
// 			  c += G_k[b1+b3*matrix_dim] * G_k_accent[b3+b2*matrix_dim];
// 			break;
			
// 		      case PARTICLE_HOLE_CHARGE:
// 			c = G_k[b1+b2*matrix_dim] * G_k_accent[b2+b1*matrix_dim];
// 			break;
		    
// 		      case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 			//cout << G_k[b1+b1*matrix_dim] << "\t" << G_k_accent[b2+b2*matrix_dim] << "\n";
// 			c = G_k[b1+b1*matrix_dim] * G_k_accent[b2+b2*matrix_dim];
// 			break;
			
// 		      default:
// 			throw std::logic_error(__FUNCTION__);
// 		      }

// 		    full_chi_0_b_k_w__b_k_w(b1, K_ind, w_vertex_index, b2, K_ind, w_vertex_index) 
// 		      += integration_factor * c;
// 		  }
// 		}
// 	      }
// 	  }
//       }

// //     {
// //       std::ofstream output_file("./chi_0_re.txt");
// //       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// // 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// // 	  output_file << real(full_chi_0_b_k_w__b_k_w(i,j)) << "\t";
// // 	}
// // 	output_file << "\n";
// //       }
// //       output_file.close();
// //     }

// //     {
// //       std::ofstream output_file("./chi_0_im.txt");
// //       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// // 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// // 	  output_file << imag(full_chi_0_b_k_w__b_k_w(i,j)) << "\t";
// // 	}
// // 	output_file << "\n";
// //       }
// //       output_file.close();
// //     }

// //     for(int i=0; i<10; i++){
// //       for(int j=0; j<10; j++){
// // 	cout << full_chi_0_b_k_w__b_k_w(i,j) << "\t";
// //       }
// //       cout << endl;
// //     }
// //     cout << endl;

//     delete [] tmp_matrix            ;
//     delete [] Sigma_matrix          ;
//     delete [] H_LDA_matrix_k        ;
//     delete [] H_LDA_matrix_k_accent ;
//     delete [] G_k                   ;
//     delete [] G_k_accent            ;
//   }



//   template<class parameter_type, class MOMS_type>
//   double analysis<parameter_type, MOMS_type>::get_integration_factor()
//   {
//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	return -1.;
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	return -2.;
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	return 1;
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }
//   }

//   template<class parameter_type, class MOMS_type>
//   std::pair<int, std::vector<double> > analysis<parameter_type, MOMS_type>::get_k_accent(int K_ind)
//   {
//     int Q_ind             = parameters.get_q_channel();
//     std::vector<double> Q = parameters.get_q_vector();//DCA_k_cluster_type::get_elements()[Q_ind];

//     std::pair<int, std::vector<double> > k_accent;

//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	{
// 	  int index                    = DCA_k_cluster_type::add(K_ind, Q_ind);
// 	  std::vector<double> K_plus_Q = DCA_k_cluster_type::get_elements()[K_ind];

// 	  for(int i=0; i<int(K_plus_Q.size()); i++)
// 	    K_plus_Q[i] += Q[i];

// 	  k_accent.first  = index;
// 	  k_accent.second = K_plus_Q;
// 	}
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	{
// 	  int index                    = DCA_k_cluster_type::add(K_ind, Q_ind);
// 	  std::vector<double> K_plus_Q = DCA_k_cluster_type::get_elements()[K_ind];

// 	  for(int i=0; i<int(K_plus_Q.size()); i++)
// 	    K_plus_Q[i] += Q[i];

// 	  k_accent.first  = index;
// 	  k_accent.second = K_plus_Q;
// 	}
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	{
// 	  int index                   = DCA_k_cluster_type::subtract(K_ind, Q_ind);
// 	  std::vector<double> Q_min_K = Q;

// 	  for(int i=0; i<int(Q_min_K.size()); i++)
// 	    Q_min_K[i] -= DCA_k_cluster_type::get_elements()[K_ind][i];

// 	  k_accent.first  = index;
// 	  k_accent.second = Q_min_K;
// 	}
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }

//     return k_accent;
//   }

//   template<class parameter_type, class MOMS_type>
//   std::pair<int, double> analysis<parameter_type, MOMS_type>::get_w_accent(int w_index, int W)
//   {
//     std::pair<int, double> w_accent;

//     switch(parameters.get_vertex_measurement_type())
//       {
//       case PARTICLE_HOLE_MAGNETIC:
// 	w_accent.first = w_index + W;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       case PARTICLE_HOLE_CHARGE:
// 	w_accent.first = w_index + W;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	w_accent.first = W + (frequency_domain_type::get_size()-1)-w_index ;
// 	assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
// 	w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
// 	break;
	
//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }

//     return w_accent;
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::load_G4_b_k_w__b_k_w()
//   {
 //    cout << __FUNCTION__ << endl;
    
//     if(VERTEX_cluster_representation == IRREDUCIBLE)
//       throw std::logic_error(__FUNCTION__);

//     assert(G4_b_k_w__b_k_w.signature() == 6);

// //     subtract_diagonal_G4_b_b_k_k_w_w();

//     int* coor_1 = new int[G4_b_k_w__b_k_w.signature()];
//     int* coor_2 = new int[MOMS.G4_k_k_w_w.signature()];
	
//     for(int i=0; i<MOMS.G4_k_k_w_w.size(); i++)
//       {
// 	MOMS.G4_k_k_w_w.linind_2_subind(i, coor_2);
 
// 	coor_1[0] = coor_2[0];
// 	coor_1[1] = coor_2[2];
// 	coor_1[2] = corresponding_extended_index[coor_2[4]];
// 	coor_1[3] = coor_2[1];
// 	coor_1[4] = coor_2[3];
// 	coor_1[5] = corresponding_extended_index[coor_2[5]];

// 	G4_b_k_w__b_k_w(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5])
// 	  = MOMS.G4_k_k_w_w(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5]);//(parameters.get_beta()*parameters.get_beta());
//       }
    
//     MOMS.G4_k_k_w_w.print_fingerprint();
//     G4_b_k_w__b_k_w.print_fingerprint();

//     {
//       std::ofstream output_file("./G4_re.txt");
//       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// 	  output_file << real(G4_b_k_w__b_k_w(i,j)) << "\t";
// 	}
// 	output_file << "\n";
//       }
//       output_file.close();
//     }

//     throw std::logic_error(__FUNCTION__);

//     {
//       std::ofstream output_file("./G4_im.txt");
//       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// 	  output_file << imag(G4_b_k_w__b_k_w(i,j)) << "\t";
// 	}
// 	output_file << "\n";
//       }
//       output_file.close();
//     }

//     interpolate_G4_b_k_w__b_k_w();

//     {
//       std::ofstream output_file("./G4_re_interp.txt");
//       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// 	  output_file << real(G4_b_k_w__b_k_w(i,j)) << "\t";
// 	}
// 	output_file << "\n";
//       }
//       output_file.close();
//     }

//     add_diagonal_G4_b_b_k_k_w_w();

//     throw std::logic_error(__FUNCTION__);

//     delete [] coor_1;
//     delete [] coor_2;
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::subtract_diagonal_G4_b_b_k_k_w_w()
//   {
//     cout << __FUNCTION__ << endl;

//     assert(b::dmn_size()==1);

//     int W        = parameters.get_number_of_positive_frequencies();
//     int W_vertex = w_VERTEX_EXTENDED::dmn_size()/2;//parameters.get_number_of_positive_vertex_frequencies();
//     int q        = parameters.get_q_channel();
    
//     int* coor = new int[3];
    
//     for(int i = 0; i < b_k_DCA_w_VERTEX_domain.get_size(); i++) 
//       {
// 	b_k_DCA_w_VERTEX_domain.linind_2_subind(i, coor);
	
// 	int b = coor[0];
// 	int k = coor[1];
// 	int w_vertex = coor[2];
// 	int w        = (corresponding_extended_index[coor[2]]- W_vertex) + W;
	    
// 	int k_plus_q  = DCA_k_cluster_type::add(k, q);     
// // 	    int q_minus_k = DCA_k_cluster_type::subtract(k, q); 

// 	switch(parameters.get_vertex_measurement_type()) 
// 	  {
// 	  case PARTICLE_HOLE_MAGNETIC:
// 	    MOMS.G4_k_k_w_w(b, b, k, k, w_vertex, w_vertex) -= -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    cout << real(MOMS.G4_k_k_w_w(b, b, k, k, w_vertex, w_vertex)) << "\t" << real(-MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w)) << endl;
// 	    //G4_0_b_k_w__b_k_w(i,i) = -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    break;
	    
// // 		  case PARTICLE_HOLE_CHARGE:
// // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = -2.*MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w);
// // 		    //G4_0_b_k_w__b_k_w(i,i) = -2.*MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// // 		    break;
		    
// // 		  case PARTICLE_PARTICLE_SUPERCONDUCTING:
// // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = MOMS.G_k_w(b1, e_UP, b1, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b2, e_UP, b2, e_UP, k, w);
// // 		    //G4_0_b_k_w__b_k_w(i,i) = MOMS.G_k_w(b, e_UP, b, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b, e_UP, b, e_UP, k, w);
// // 		    //cout << frequency_domain_type::get_elements()[2*W-1-w] << "  <==>  " << frequency_domain_type::get_elements()[w] << endl;
// // 		    break;
		    
// 	  default:
// 	    throw std::logic_error(__FUNCTION__);
// 	  }
//       }
  
// //     throw std::logic_error("STOP");

//     delete [] coor;
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::add_diagonal_G4_b_b_k_k_w_w()
//   {
//     assert(b::dmn_size()==1);

//     int W        = parameters.get_number_of_positive_frequencies();
//     int W_vertex = w_VERTEX_EXTENDED::dmn_size()/2;//parameters.get_number_of_positive_vertex_frequencies();
//     int q        = parameters.get_q_channel();
    
//     int* coor = new int[3];
    
//     for(int i = 0; i < b_k_DCA_w_VERTEX_domain.get_size(); i++) 
//       {
// 	b_k_DCA_w_VERTEX_domain.linind_2_subind(i, coor);
	
// 	int b = coor[0];
// 	int k = coor[1];
// 	int w_vertex = coor[2];
// 	int w        = (corresponding_extended_index[coor[2]]- W_vertex) + W;
	    
// 	int k_plus_q  = DCA_k_cluster_type::add(k, q);     
// // 	    int q_minus_k = DCA_k_cluster_type::subtract(k, q); 

// 	switch(parameters.get_vertex_measurement_type()) 
// 	  {
// 	  case PARTICLE_HOLE_MAGNETIC:
// 	    MOMS.G4_k_k_w_w(b, b, k, k, w_vertex, w_vertex) += -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    //cout << real(G4_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex)) << "\t" << real(-MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w)) << endl;
// 	    //G4_0_b_k_w__b_k_w(i,i) = -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    break;
	    
// // 		  case PARTICLE_HOLE_CHARGE:
// // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = -2.*MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w);
// // 		    //G4_0_b_k_w__b_k_w(i,i) = -2.*MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// // 		    break;
		    
// // 		  case PARTICLE_PARTICLE_SUPERCONDUCTING:
// // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = MOMS.G_k_w(b1, e_UP, b1, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b2, e_UP, b2, e_UP, k, w);
// // 		    //G4_0_b_k_w__b_k_w(i,i) = MOMS.G_k_w(b, e_UP, b, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b, e_UP, b, e_UP, k, w);
// // 		    //cout << frequency_domain_type::get_elements()[2*W-1-w] << "  <==>  " << frequency_domain_type::get_elements()[w] << endl;
// // 		    break;
		    
// 	  default:
// 	    throw std::logic_error(__FUNCTION__);
// 	  }
//       }
  
// //     throw std::logic_error("STOP");

//     delete [] coor;
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::interpolate_G4_b_k_w__b_k_w()
//   {
//     cout << __FUNCTION__ << endl;

// //     subtract_diagonal_G4_b_k_w__b_k_w();

//     int* coor_1 = new int[G4_b_k_w__b_k_w.signature()];
// //     int* coor_2 = new int[MOMS.G4_k_k_w_w.signature()];
	
//     for(int i=0; i<G4_b_k_w__b_k_w.size(); i++)
//       {
// 	G4_b_k_w__b_k_w.linind_2_subind(i, coor_1);
 
// 	double w1 = w_VERTEX_EXTENDED::parameter_type::get_elements()[coor_1[2]];
//  	double w2 = w_VERTEX_EXTENDED::parameter_type::get_elements()[coor_1[5]];

// 	int w1_lb = 0;
// 	while(w_VERTEX::parameter_type::get_elements()[w1_lb+1] < w1 && w1_lb<w_VERTEX::dmn_size()-2)
// 	  w1_lb++;

// 	int w1_ub = w1_lb+1;

// 	int w2_lb = 0;
// 	while(w_VERTEX::parameter_type::get_elements()[w2_lb+1] < w2 && w2_lb<w_VERTEX::dmn_size()-2)
// 	  w2_lb++;

// 	int w2_ub = w2_lb+1;

// // 	cout << w_VERTEX::parameter_type::get_elements()[w1_lb] << "\t" <<  w1 << "\t" <<  w_VERTEX::parameter_type::get_elements()[w1_ub] << "\n";
// //  	cout << w_VERTEX::parameter_type::get_elements()[w2_lb] << "\t" <<  w2 << "\t" <<  w_VERTEX::parameter_type::get_elements()[w2_ub] << "\n";
// // 	cout << "\n";

// 	assert(w_VERTEX::parameter_type::get_elements()[w1_lb] <= w1 && w1 <= w_VERTEX::parameter_type::get_elements()[w1_ub]);
//  	assert(w_VERTEX::parameter_type::get_elements()[w2_lb] <= w2 && w2 <= w_VERTEX::parameter_type::get_elements()[w2_ub]);

// 	if(abs(G4_b_k_w__b_k_w(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5])) < 1.e-12)
// 	  {
// 	    std::complex<double> f_0_0 = MOMS.G4_k_k_w_w(coor_1[0], coor_1[3], coor_1[1], coor_1[4], w1_lb , w2_lb);//(parameters.get_beta()*parameters.get_beta()); 
// 	    std::complex<double> f_1_0 = MOMS.G4_k_k_w_w(coor_1[0], coor_1[3], coor_1[1], coor_1[4], w1_ub , w2_lb);//(parameters.get_beta()*parameters.get_beta()); 
// 	    std::complex<double> f_0_1 = MOMS.G4_k_k_w_w(coor_1[0], coor_1[3], coor_1[1], coor_1[4], w1_lb , w2_ub);//(parameters.get_beta()*parameters.get_beta()); 
// 	    std::complex<double> f_1_1 = MOMS.G4_k_k_w_w(coor_1[0], coor_1[3], coor_1[1], coor_1[4], w1_ub , w2_ub);//(parameters.get_beta()*parameters.get_beta()); 

// 	    double x_0 = w_VERTEX::parameter_type::get_elements()[w1_lb];
// 	    double x_1 = w_VERTEX::parameter_type::get_elements()[w1_ub];
// 	    double y_0 = w_VERTEX::parameter_type::get_elements()[w2_lb];
// 	    double y_1 = w_VERTEX::parameter_type::get_elements()[w2_ub];
	    
// 	    double t0 = (w1-x_0)/(x_1-x_0);
// 	    double t1 = (w2-y_0)/(y_1-y_0);

// 	    G4_b_k_w__b_k_w(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5]) = (1-t0)*(1-t1)*f_0_0 + t0*(1-t1)*f_1_0 + (1-t0)*t1*f_0_1 + t0*t1*f_1_1;
// 	  }
//       }

//     add_diagonal_G4_b_k_w__b_k_w();
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::add_diagonal_G4_b_k_w__b_k_w()
//   {
//     int W        = parameters.get_number_of_positive_frequencies();
//     int W_vertex = w_VERTEX_EXTENDED::dmn_size()/2;//parameters.get_number_of_positive_vertex_frequencies();
//     int q        = parameters.get_q_channel();
    
//     int* coor = new int[3];
    
//     for(int i = 0; i < b_k_DCA_w_VERTEX_EXTENDED_domain.get_size(); i++) 
//       {
// 	b_k_DCA_w_VERTEX_EXTENDED_domain.linind_2_subind(i, coor);
	
// 	int b = coor[0];
// 	int k        = coor[1];
// 	int w_vertex = coor[2];
// 	int w        = (coor[2] - W_vertex) + W;
	    
// 	int k_plus_q  = DCA_k_cluster_type::add(k, q);     
// 	// 	    int q_minus_k = DCA_k_cluster_type::subtract(k, q); 
	
// 	switch(parameters.get_vertex_measurement_type()) 
// 	  {
// 	  case PARTICLE_HOLE_MAGNETIC:
	    
// 	    G4_b_k_w__b_k_w(b,k,w_vertex,b,k,w_vertex) += -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    //cout << real(G4_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex)) << "\t" << real(-MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w)) << endl;
// 	    //G4_0_b_k_w__b_k_w(i,i) = -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    break;
	    
// 	    // 		  case PARTICLE_HOLE_CHARGE:
// 	    // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = -2.*MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w);
// 	    // 		    //G4_0_b_k_w__b_k_w(i,i) = -2.*MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 	    // 		    break;
	    
// 	    // 		  case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 	    // 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = MOMS.G_k_w(b1, e_UP, b1, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b2, e_UP, b2, e_UP, k, w);
// 	    // 		    //G4_0_b_k_w__b_k_w(i,i) = MOMS.G_k_w(b, e_UP, b, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b, e_UP, b, e_UP, k, w);
// 	    // 		    //cout << frequency_domain_type::get_elements()[2*W-1-w] << "  <==>  " << frequency_domain_type::get_elements()[w] << endl;
// 	    // 		    break;
	    
// 	  default:
// 	    throw std::logic_error(__FUNCTION__);
// 	  }
//       }
//     //     throw std::logic_error("STOP");

//     delete [] coor;
//   }


//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::load_G4_b_k_w__b_k_w()
//   {
//     cout << __FUNCTION__ << endl;

//     assert(G4_b_k_w__b_k_w.signature() == 6);

//     int* coor_1 = new int[G4_b_k_w__b_k_w.signature()];
//     int* coor_2 = new int[MOMS.G4_k_k_w_w.signature()];
	
//     for(int i=0; i<G4_b_k_w__b_k_w.size(); i++)
//       {
// 	G4_b_k_w__b_k_w.linind_2_subind(i, coor_1);
 
// 	coor_2[0] = coor_1[0];
// 	coor_2[1] = coor_1[3];
// 	coor_2[2] = coor_1[1];
// 	coor_2[3] = coor_1[4];
// 	coor_2[4] = coor_1[2];
// 	coor_2[5] = coor_1[5];

// 	if(VERTEX_cluster_representation == IRREDUCIBLE)
// 	    throw std::logic_error(__FUNCTION__);

// 	G4_b_k_w__b_k_w(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5])
// 	  = MOMS.G4_k_k_w_w(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5])/(parameters.get_beta()*parameters.get_beta());
//       }
    
// //     {
// //       std::ofstream output_file("./G4_re.txt");
// //       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// // 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// // 	  output_file << real(G4_b_k_w__b_k_w(i,j)) << "\t";
// // 	}
// // 	output_file << "\n";
// //       }
// //       output_file.close();
// //     }

// //     {
// //       std::ofstream output_file("./G4_im.txt");
// //       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// // 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// // 	  output_file << imag(G4_b_k_w__b_k_w(i,j)) << "\t";
// // 	}
// // 	output_file << "\n";
// //       }
// //       output_file.close();
// //     }

//     delete [] coor_1;
//     delete [] coor_2;
//   }

//   template<class parameter_type, class MOMS_type>
//   void analysis<parameter_type, MOMS_type>::load_G4_0_b_k_w__b_k_w()
//   {
//     cout << __FUNCTION__ << endl;

//     int W        = parameters.get_number_of_positive_frequencies();
//     int W_vertex = w_VERTEX_EXTENDED::dmn_size()/2;//parameters.get_number_of_positive_vertex_frequencies();
//     int q        = parameters.get_q_channel();

//     int* coor = new int[3];

//     for(int i = 0; i < b_k_DCA_w_VERTEX_EXTENDED_domain.get_size(); i++) 
//       {
// 	b_k_DCA_w_VERTEX_EXTENDED_domain.linind_2_subind(i, coor);
	
// 	int b = coor[0];

// 	if(b == 0)
// 	  {
// 	    int k        = coor[1];
// 	    int w_vertex = coor[2];
// 	    int w = (coor[2] - W_vertex) + W;
	    
// 	    int k_plus_q  = DCA_k_cluster_type::add(k, q);     
// 	    int q_minus_k = DCA_k_cluster_type::subtract(k, q); 
	    
// 	    for(int b1=0; b1<b::dmn_size(); b1++){
// 	      for(int b2=0; b2<b::dmn_size(); b2++){

// 		switch(parameters.get_vertex_measurement_type()) 
// 		  {
// 		  case PARTICLE_HOLE_MAGNETIC:
// 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = -MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w);
// 		    //G4_0_b_k_w__b_k_w(i,i) = -MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 		    break;
		
// 		  case PARTICLE_HOLE_CHARGE:
// 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = -2.*MOMS.G_k_w(b1, e_UP, b2, e_UP, k, w)*MOMS.G_k_w(b2, e_UP, b1, e_UP, k_plus_q, w);
// 		    //G4_0_b_k_w__b_k_w(i,i) = -2.*MOMS.G_k_w(b, e_UP, b, e_UP, k, w)*MOMS.G_k_w(b, e_UP, b, e_UP, k_plus_q, w);
// 		    break;
		    
// 		  case PARTICLE_PARTICLE_SUPERCONDUCTING:
// 		    G4_0_b_k_w__b_k_w(b1,k,w_vertex,b2,k,w_vertex) = MOMS.G_k_w(b1, e_UP, b1, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b2, e_UP, b2, e_UP, k, w);
// 		    //G4_0_b_k_w__b_k_w(i,i) = MOMS.G_k_w(b, e_UP, b, e_UP, q_minus_k, 2*W-1-w)*MOMS.G_k_w(b, e_UP, b, e_UP, k, w);
// 		    //cout << frequency_domain_type::get_elements()[2*W-1-w] << "  <==>  " << frequency_domain_type::get_elements()[w] << endl;
// 		    break;
		    
// 		  default:
// 		    throw std::logic_error(__FUNCTION__);
// 		  }
// 	      }
// 	    }
// 	  }
//       }

//     delete [] coor;

//     {
//       std::ofstream output_file("./G4_0_re.txt");
//       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// 	  output_file << real(G4_0_b_k_w__b_k_w(i,j)) << "\t";
// 	}
// 	output_file << "\n";
//       }
//       output_file.close();
//     }

//     {
//       std::ofstream output_file("./G4_0_im.txt");
//       for(int i=0; i<std::sqrt(G4_b_k_w__b_k_w.size()); i++){
// 	for(int j=0; j<std::sqrt(G4_b_k_w__b_k_w.size()); j++){
// 	  output_file << imag(G4_0_b_k_w__b_k_w(i,j)) << "\t";
// 	}
// 	output_file << "\n";
//       }
//       output_file.close();
//     }

//    }


/*
  template<class parameter_type, class MOMS_type>
  std::complex<double> analysis<parameter_type, MOMS_type>::compute_full_chi_0_tail()
  {
    cout << __FUNCTION__ << endl;

    std::complex<double> result = 0.;

    int                 W     = parameters.get_w_channel();
    double chemical_potential = parameters.get_chemical_potential();

    const static double               Nb_interpolation = 16;
    std::vector<std::vector<double> > centered_mesh    = Mesh<DCA_k_cluster_type>::execute(int(Nb_interpolation));
    std::vector<std::vector<double> > mesh             = centered_mesh;

    double integration_factor = get_integration_factor()/double(centered_mesh.size());

    int matrix_size = 2*2*b::dmn_size()*b::dmn_size();
    int matrix_dim  = 2*b::dmn_size();
    
    std::complex<double>* tmp_matrix            = new std::complex<double>[matrix_size];
    std::complex<double>* Sigma_matrix          = new std::complex<double>[matrix_size];
    std::complex<double>* H_LDA_matrix_k        = new std::complex<double>[matrix_size];
    std::complex<double>* H_LDA_matrix_k_accent = new std::complex<double>[matrix_size];
    std::complex<double>* G_k                   = new std::complex<double>[matrix_size];
    std::complex<double>* G_k_accent            = new std::complex<double>[matrix_size];
    
    wannier_interpolation::mesh_k::get_size() = int(mesh.size());
    wannier_interpolation WIP_k       (MOMS.H_LDA);
    wannier_interpolation WIP_k_accent(MOMS.H_LDA);

    invert_plan<std::complex<double> > invert_pln(matrix_dim);

    for(int K_ind=0; K_ind<DCA_k_cluster_type::get_size(); K_ind++)
      {
	std::vector<double> K = DCA_k_cluster_type::get_elements()[K_ind];
	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K);
	wannier_interpolation::H_k_interpolated_type& H_k        = WIP_k.execute(mesh);

	int                 K_accent_ind = get_k_accent(K_ind).first;
	std::vector<double> K_accent     = get_k_accent(K_ind).second;
	Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K_accent);
	wannier_interpolation::H_k_interpolated_type& H_k_accent = WIP_k_accent.execute(mesh);

	// \sum_{k \in K_{k}}
	for(int k_ind=0; k_ind<int(mesh.size()); k_ind++)
	  {
	    memcpy(H_LDA_matrix_k       , &H_k       (0,0,k_ind), sizeof(std::complex<double>)*matrix_size);
	    memcpy(H_LDA_matrix_k_accent, &H_k_accent(0,0,k_ind), sizeof(std::complex<double>)*matrix_size);

	    for(int w_index=w_VERTEX::dmn_size(); w_index<w::dmn_size()/2; w_index++)
	      {
		double w  = w::parameter_type::get_elements()[w_index+ w::dmn_size()/2];
		int w_ind = w_index  + w::dmn_size()/2;
		
		{// G(k)
		  memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_ind, w_ind), sizeof(std::complex<double>)*matrix_size);

		  {// G^-1 =  -(H_k + Sigma) + i*w + mu
		    for(int index=0; index < matrix_size; index++)
		      tmp_matrix[index] = -(H_LDA_matrix_k[index] + Sigma_matrix[index]);  
		    
		    for(int nu=0; nu<matrix_dim; nu++)
		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w);
		  }
		  
		  {// G(k)^-1 = (-(H(k)+Sigma(k)) + i*w) ==> G(k) = (-(H(k)+Sigma(k)) + i*w)^-1 
		    memcpy(invert_pln.Matrix, tmp_matrix                    , sizeof(std::complex<double>)*matrix_size);
		    invert_pln.execute_plan();
		    memcpy(&G_k[0]          , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
		  }
		}

		{// G(k+Q)
		  std::pair<int, double> w_accent = get_w_accent(w_ind, W);
		  
		  memcpy(Sigma_matrix, &MOMS.Sigma(0,0,K_accent_ind, w_accent.first), sizeof(std::complex<double>)*matrix_size);

		  {// G^-1 =  -(H_k_accent + Sigma) + i*w + mu
		    for(int index=0; index < matrix_size; index++)
		      tmp_matrix[index] = -(H_LDA_matrix_k_accent[index] + Sigma_matrix[index]);  
		    
		    for(int nu=0; nu<matrix_dim; nu++)
		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w_accent.second);
		  }
		  
		  {// G(k+Q)^-1 = (-(H(k+Q)+Sigma(k+Q)) + i*w) ==> G(k+Q) = (-(H(k+Q)+Sigma(k+Q)) + i*w)^-1 
		    memcpy(invert_pln.Matrix, tmp_matrix                     , sizeof(std::complex<double>)*matrix_size);
		    invert_pln.execute_plan();
		    memcpy(&G_k_accent[0]    , &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
		  }
		}

		for(int b1=0; b1<b::dmn_size(); b1++){
		  for(int b2=0; b2<b::dmn_size(); b2++){

		    std::complex<double> c(0,0);

		    switch(parameters.get_vertex_measurement_type()) 
		      {
		      case PARTICLE_HOLE_MAGNETIC:
			for(int b3=0; b3<b::dmn_size(); b3++)
			  c += G_k[b1+b3*matrix_dim] * G_k_accent[b3+b2*matrix_dim];
			break;
			
		      case PARTICLE_HOLE_CHARGE:
			c = G_k[b1+b2*matrix_dim] * G_k_accent[b2+b1*matrix_dim];
			break;
		    
		      case PARTICLE_PARTICLE_SUPERCONDUCTING:
			//cout << G_k[b1+b1*matrix_dim] << "\t" << G_k_accent[b2+b2*matrix_dim] << "\n";
			c = G_k[b1+b1*matrix_dim] * G_k_accent[b2+b2*matrix_dim];
			break;
			
		      default:
			throw std::logic_error(__FUNCTION__);
		      }

		    result += integration_factor * (c + std::conj(c));
		  }
		}
	      }
	  }
      }

    delete [] tmp_matrix            ;
    delete [] Sigma_matrix          ;
    delete [] H_LDA_matrix_k        ;
    delete [] H_LDA_matrix_k_accent ;
    delete [] G_k                   ;
    delete [] G_k_accent            ;

    return result;
  }
 */

  /*    
  template<class parameter_type, class MOMS_type>
  void analysis<parameter_type, MOMS_type>::compute_full_chi_0_b_k_w__b_k_w()
  {
    cout << __FUNCTION__ << endl;

    int                 W = parameters.get_w_channel();
    std::vector<double> Q = DCA_k_cluster_type::get_elements()[parameters.get_q_channel()];
    double chemical_potential = parameters.get_chemical_potential();

    const static double               Nb_interpolation = 100.;
    std::vector<std::vector<double> > centered_mesh    = Mesh<DCA_k_cluster_type>::execute(Nb_interpolation);
    std::vector<std::vector<double> > mesh             = centered_mesh;

    double integration_factor = get_integration_factor()/double(centered_mesh.size());

    //const static double Nb_interpolation = 10.;

    int matrix_size = 2*2*b::dmn_size()*b::dmn_size();
    int matrix_dim  = 2*b::dmn_size();
    
    static std::complex<double>* tmp_matrix   = new std::complex<double>[matrix_size];
    static std::complex<double>* Sigma_matrix = new std::complex<double>[matrix_size];
    static std::complex<double>* H_LDA_matrix = new std::complex<double>[matrix_size];
    static std::complex<double>* G_k          = new std::complex<double>[matrix_size];
    static std::complex<double>* G_k_accent   = new std::complex<double>[matrix_size];
    static int*                 coor          = new int[6];
    memset(coor, 0, sizeof(int)*6);
    
    invert_plan<std::complex<double> > invert_pln(matrix_dim);

    for(int i=0; i<w_VERTEX::dmn_size(); i++)
      {
	double w = w_VERTEX::parameter_type::get_elements()[i];
	coor[5] = i - w_VERTEX::dmn_size()/2 + w::dmn_size()/2;

	assert( fabs(w -  w::parameter_type::get_elements()[coor[5]]) < 1.e-6);

	for(int j=0; j<DCA_k_cluster_type::get_size(); j++)
	  {
	    coor[4] = j;

	    std::vector<double> K = DCA_k_cluster_type::get_elements()[j];
	    Mesh<DCA_k_cluster_type>::translate_mesh(centered_mesh, mesh, K);

	    for(int mesh_ind=0; mesh_ind<int(mesh.size()); mesh_ind++)
	      {
		{
		  //cout << w << "\t";

		  memset(coor, 0, sizeof(int)*4);
		  interpolation<LDA_k_cluster_type, 1>::execute_on_matrix(&mesh[mesh_ind][0], MOMS.H_LDA, coor, matrix_size, H_LDA_matrix);
		  
		  memset(coor, 0, sizeof(int)*4);
		  interpolation<DCA_k_cluster_type, 0>::execute_on_matrix(&mesh[mesh_ind][0], MOMS.Sigma, coor, matrix_size, Sigma_matrix);
		  
		  {// G^-1 =  -(H + Sigma) + i*w
		    for(int index=0; index < matrix_size; index++)
		      tmp_matrix[index] = -(H_LDA_matrix[index] + Sigma_matrix[index]);  
		    
		    for(int nu=0; nu<matrix_dim; nu++)
		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w);
		  }
		  
		  {// G(k)^-1 = (-(H(k)+Sigma(k)) + i*w) ==> G(k) = (-(H(k)+Sigma(k)) + i*w)^-1 
		    invert_pln.Matrix = tmp_matrix;
		    invert_pln.execute_plan();
		    memcpy(&G_k[0], &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
		  }
		}

		{
		  std::vector<double>    k_accent = get_k_accent(mesh[mesh_ind], Q);
		  std::pair<int, double> w_accent = get_w_accent(coor[5]       , W);

		  //cout << w_accent.second << endl;

		  memset(coor, 0, sizeof(int)*4);
		  interpolation<LDA_k_cluster_type, 1>::execute_on_matrix(&k_accent[0], MOMS.H_LDA, coor, matrix_size, H_LDA_matrix);
		  
		  memset(coor, 0, sizeof(int)*4);
		  interpolation<DCA_k_cluster_type, 0>::execute_on_matrix(&k_accent[0], MOMS.Sigma, coor, matrix_size, Sigma_matrix);
		  		
		  {// G^-1 =  -(H + Sigma) + i*w
		    for(int index=0; index < matrix_size; index++)
		      tmp_matrix[index] = -(H_LDA_matrix[index] + Sigma_matrix[index]);  
		    
		    for(int nu=0; nu<matrix_dim; nu++)
		      tmp_matrix[nu + matrix_dim*nu] += std::complex<double>(chemical_potential, w_accent.second);
		  }
		  
		  {// G(k+Q)^-1 = (-(H(k+Q)+Sigma(k+Q)) + i*w) ==> G(k+Q) = (-(H(k+Q)+Sigma(k+Q)) + i*w)^-1 
		    invert_pln.Matrix = tmp_matrix;
		    invert_pln.execute_plan();
		    memcpy(&G_k_accent[0], &invert_pln.inverted_matrix[0], sizeof(std::complex<double>)*matrix_size);
		  }
		}

		for(int b1=0; b1<b::dmn_size(); b1++){
		  for(int b2=0; b2<b::dmn_size(); b2++){

		    std::complex<double> c;

		    switch(parameters.get_vertex_measurement_type()) 
		      {
		      case PARTICLE_HOLE_MAGNETIC:
			c = G_k[b1+b2*matrix_dim] * G_k_accent[b2+b1*matrix_dim];
			break;
			
		      case PARTICLE_HOLE_CHARGE:
			c = G_k[b1+b2*matrix_dim] * G_k_accent[b2+b1*matrix_dim];
			break;
		    
		      case PARTICLE_PARTICLE_SUPERCONDUCTING:
			c = G_k[b1+b1*matrix_dim] * G_k_accent[b2+b2*matrix_dim];
			break;
			
		      default:
			throw std::logic_error(__FUNCTION__);
		      }

		    full_chi_0_b_k_w__b_k_w(b1, j, i, b2, j, i) 
		      += integration_factor * c;//G_k[b1+b2*matrix_dim] * G_k_accent[b2+b1*matrix_dim];

		  }
		}
	      }
	  }
      }

//     int* coor_1 = new int[G4_b_k_w__b_k_w.signature()];
//     int* coor_2 = new int[MOMS.G4_k_k_w_w.signature()];
    
//     for(int i=0; i<10; i++){
//       full_chi_0_b_k_w__b_k_w.linind_2_subind(i, coor_1);
//       for(int j=0; j<10; j++){
// 	full_chi_0_b_k_w__b_k_w.linind_2_subind(j, coor_2);
// 	//if(coor_1[0] == 0 && coor_2[0] == 0)
// 	  cout << real(full_chi_0_b_k_w__b_k_w(i,j)) << "\t";
//       }
//       cout << endl;
//     }
//     cout << endl;
    
//     for(int i=0; i<10; i++){
//       full_chi_0_b_k_w__b_k_w.linind_2_subind(i, coor_1);
//       for(int j=0; j<10; j++){
// 	full_chi_0_b_k_w__b_k_w.linind_2_subind(j, coor_2);
// 	if(coor_1[0] == 1 && coor_2[0] == 1)
// 	  cout << real(full_chi_0_b_k_w__b_k_w(i,j)) << "\t";
//       }
//       cout << endl;
//     }
//     cout << endl;

    //delete [] coor_1;
    //delete [] coor_2;

  }
*/

