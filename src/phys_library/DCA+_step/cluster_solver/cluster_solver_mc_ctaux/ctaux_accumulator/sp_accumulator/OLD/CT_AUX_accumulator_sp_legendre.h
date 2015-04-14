//-*-C++-*-

#ifndef CT_AUX_ACCUMULATOR_SP_LEGENDRE_H
#define CT_AUX_ACCUMULATOR_SP_LEGENDRE_H

namespace QMC {
  
  /*! 
   *  \ingroup CT-AUX
   *
   *  \brief   This class measures the single-particle functions with an LEGENDRE scheme in the CT-AUX QMC engine.
   *  \author  Peter Staar
   *  \version 1.0
   */
  template<class parameters_type>
  class MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>
  {
#include "type_definitions.h" 

    typedef double scalar_type;

    typedef vertex_singleton<base_cluster_type> vertex_singleton_type;

//     typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
//     typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

//     typedef dmn_0<r_cluster_type> r_dmn_t;
//     typedef dmn_0<k_cluster_type> k_dmn_t;

    typedef r_DCA r_dmn_t;
    typedef k_DCA k_dmn_t;

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::Concurrency_Type concurrency_type;

//     typedef legendre_domain<time_domain, SINGLE_PARTICLE_QUANTITY>  legendre_t;
//     typedef dmn_0<legendre_t>                                       legendre_sp_dmn_t;

  public:

    MC_single_particle_accumulator(parameters_type& parameters_ref);

    ~MC_single_particle_accumulator();

    void initialize(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w);
    
    void finalize();

    /*!
     *  \brief Print the functions M_r_w and M_k_w.
     */
    template<class stream_type>
    void to_JSON(stream_type& ss);

    template<class configuration_type, class vertex_vertex_matrix_type>
    void accumulate_M_r_w(configuration_type&        configuration_e_spin,
			  vertex_vertex_matrix_type& M,
			  double                     sign,
			  e_spin_states              e_spin);

    void compute_M_r_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    std::vector<double> P;

    FUNC_LIB::function<double, dmn_4<legendre_sp_dmn_t,nu,nu,r_dmn_t> > M_l_r;
    FUNC_LIB::function<double, dmn_4<nu,nu,r_dmn_t,legendre_sp_dmn_t> > M_r_l;

    FUNC_LIB::function<double, dmn_4<nu,nu,r_dmn_t,t> > M_r_t;
  };

  template<class parameters_type>
  MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::MC_single_particle_accumulator(parameters_type& parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency()),

    P(parameters.get_nb_of_legendre_coefficients_single_particle(), 0),

    M_l_r("M_legendre_r"),
    M_r_l("M_r_legendre"),

    M_r_t("M_r_t")
  {}

  template<class parameters_type>
  MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::~MC_single_particle_accumulator()
  {}

  template<class parameters_type>
  void MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::initialize(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>  >& M_r_w)
  {
    for(int i=0; i<M_l_r.size(); i++)
      M_l_r(i) = 0;

    for(int i=0; i<M_r_l.size(); i++)
      M_r_l(i) = 0;

    for(int i=0; i<M_r_w.size(); i++)
      M_r_w(i) = 0;
  }

  template<class parameters_type>
  void MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::finalize()
  {
    concurrency.sum_and_average(M_r_l, parameters.get_number_of_measurements());

    if(concurrency.id()==0){
      for(int l=0; l<legendre_sp_dmn_t::dmn_size(); l++){
	cout << "\t" << l;
	for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
	  cout << "\t" << M_r_l(0,0,r_ind,l);
	}
	cout << "\n";
      }
    }
  }

  template<class parameters_type>
  template<class stream_type>
  void MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::to_JSON(stream_type& ss)
  {
    ss << ",";
    M_r_l.to_JSON(ss);

    ss << ",";
    M_r_t.to_JSON(ss);
  }

  template<class parameters_type>
  template<class configuration_type, class vertex_vertex_matrix_type>
  void MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::accumulate_M_r_w(configuration_type&        configuration_e_spin,
											    vertex_vertex_matrix_type& M,
											    double                     sign,
											    e_spin_states              e_spin)
  {
    scalar_type  beta         = parameters.get_beta();
    scalar_type one_div_beta = 1./parameters.get_beta();
    
    int spin_index = do_cast<int>::execute(e_spin);

    int         r_ind, b_i, b_j, r_i, r_j, s_i, s_j;
    scalar_type t_i, t_j, delta_tau, scaled_tau, f_tau;

    int configuration_size = configuration_e_spin.size();

    for(int j=0; j<configuration_size; j++)
      {
	vertex_singleton_type& configuration_e_spin_j = configuration_e_spin[j];

	b_j = configuration_e_spin_j.get_band();
	r_j = configuration_e_spin_j.get_r_site();
	t_j = configuration_e_spin_j.get_tau();
	s_j = configuration_e_spin_j.get_HS_spin();

	for(int i=0; i<configuration_size; i++)
	  {
	    if(configuration_e_spin[i].get_HS_spin() != HS_ZERO && configuration_e_spin[j].get_HS_spin() != HS_ZERO)
	      {
		vertex_singleton_type& configuration_e_spin_i = configuration_e_spin[i];
		
		b_i = configuration_e_spin_i.get_band();
		r_i = configuration_e_spin_i.get_r_site();
		t_i = configuration_e_spin_i.get_tau();
		s_i = configuration_e_spin_i.get_HS_spin();
		
		r_ind = r_cluster_type::subtract(r_j, r_i);
		
		delta_tau = t_i-t_j;
		
		if(delta_tau<0){
		  scaled_tau = 2*(delta_tau+beta)*one_div_beta-1.;
		  f_tau      = -M(i,j)*sign;
		}
		else{
		  scaled_tau = 2*(delta_tau)*one_div_beta-1.;
		  f_tau      = M(i,j)*sign;
		}

		assert(configuration_e_spin[i].get_HS_spin() != HS_ZERO 
		       && configuration_e_spin[j].get_HS_spin() != HS_ZERO
		       && square(s_i*s_j) == 1);
		
		legendre_basis_functions<time_domain>::evaluate_more_than_2(scaled_tau, P);
		
		scalar_type* ptr = &M_l_r(0, b_i, spin_index, b_j, spin_index, r_ind);

		for(int l=0; l<legendre_sp_dmn_t::dmn_size(); l++)
		  ptr[l] += f_tau*P[l];
	      }
	  }
      }
  }      
  
  template<class parameters_type>
  void MC_single_particle_accumulator<CT_AUX, LEGENDRE, parameters_type>::compute_M_r_w(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w)
  {
    cout << scientific;
    cout.precision(6);

    for(int l=0; l<legendre_sp_dmn_t::dmn_size(); l++)
      for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++)
	for(int nu_j=0; nu_j<b::dmn_size()*s::dmn_size(); nu_j++)
	  for(int nu_i=0; nu_i<b::dmn_size()*s::dmn_size(); nu_i++)
	    M_r_l(nu_i,nu_j,r_ind,l) = M_l_r(l,nu_i,nu_j,r_ind);

    legendre_transformation<legendre_sp_dmn_t, t>::execute(M_r_l, M_r_t);

    for(int l=0; l<legendre_sp_dmn_t::dmn_size(); l++){

      cout << l;

      double renorm = (2*l+1)/parameters.get_beta();//sqrt(2*l+1)/parameters.get_beta();

      for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
	for(int nu_j=0; nu_j<b::dmn_size()*s::dmn_size(); nu_j++){
	  for(int nu_i=0; nu_i<b::dmn_size()*s::dmn_size(); nu_i++){
	    M_r_l(nu_i,nu_j,r_ind,l) = renorm*M_l_r(l,nu_i,nu_j,r_ind);
	  }
	}

	cout << "\t" << M_r_l(0,0,r_ind,l);
      }
      cout << "\n";
    }

    //legendre_transformation<legendre_sp_dmn_t, t>::execute(M_r_l, M_r_t);

    FT<legendre_sp_dmn_t, w>::execute(M_r_l, M_r_w);

    double one_div_n_sites = 1./double(base_cluster_type::get_cluster_size());
    M_r_w *= one_div_n_sites;
  }

}
#endif
