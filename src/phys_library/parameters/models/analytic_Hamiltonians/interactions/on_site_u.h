//-*-C++-*-

#ifndef ON_SITE_U_H
#define ON_SITE_U_H

/*!
 *   \ingroup MODELS
 *
 *   \author Peter Staar
 *   \brief  Simple template that implements on-site interaction, for different MC-solvers. 
 */
class on_site_u
{
public:

  template<class vertex_pair_type, class parameters_type, class concurrency_type, class H_interactions_function_type>
  static void set_vertex(vertex_pair_type& vertex, 
			 parameters_type&  parameters, 
			 concurrency_type& concurrency,
			 H_interactions_function_type& H_interactions)
  {
    typedef typename H_interactions_function_type::this_domain_type interaction_domain_type;

    typedef typename interaction_domain_type::domain_typelist_0 dmn_typelist_0;
    typedef typename interaction_domain_type::domain_typelist_2 dmn_typelist_2;

//     typedef typename TypeAt<dmn_typelist_0, 0>::Result band_domain_type;
    typedef typename TypeAt<dmn_typelist_0, 1>::Result spin_domain_type;

    typedef typename TypeAt<dmn_typelist_2, 0>::Result r_dmn_type;

    typedef typename parameters_type::nu nu;

    int BANDS = parameters.get_interacting_bands().size();

    do
      {
	int band_ind_1 = concurrency.get_random_number()*BANDS;
	//int band_ind_2 = concurrency.get_random_number()*BANDS;

        // In order to not double-count the intra-band interaction we need different probabilities
        // for chosing the second band to be the same/different to the first band. Each band
        // different to the first band should have twice the probablility to be chosen than the same
        // band.
        int band_ind_2 = band_ind_1;

        // Efficient version for 2 bands.
        // FIXME: Change to general version later !!!
        assert(BANDS == 2);
        double p = 2./3;
        if (concurrency.get_random_number() < p) {
          band_ind_2 = 1 - band_ind_1;
        }

        // General version.
        // double p = 1. - 1./(2.*BANDS - 1.);  // probablity to choose any different band
        // if (concurrency.get_random_number() < p) {
        //   while (band_ind_2 == band_ind_1) {
        //     band_ind_2 = concurrency.get_random_number()*BANDS;
        //   }
        // }
      
	vertex.get_bands().first  = parameters.get_interacting_bands()[band_ind_1];
	vertex.get_bands().second = parameters.get_interacting_bands()[band_ind_2];
      
	vertex.get_e_spins().first = spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	if(vertex.get_bands().first == vertex.get_bands().second)
	  {
	    if(vertex.get_e_spins().first == e_UP)
	      vertex.get_e_spins().second = e_DN;
	    else
	      vertex.get_e_spins().second = e_UP;
	  }
	else
	  vertex.get_e_spins().second = spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	vertex.get_spin_orbitals().first  = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().first , vertex.get_e_spins().first); 
	vertex.get_spin_orbitals().second = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().second, vertex.get_e_spins().second); 
      }
    while(std::fabs(H_interactions(vertex.get_spin_orbitals().first, vertex.get_spin_orbitals().second, 0)) < 1.e-3 );

    int r_site = int(r_dmn_type::get_size()*concurrency.get_random_number());

//     cout << r_site << "\t" << vertex.get_spin_orbitals().first << "\t" << vertex.get_spin_orbitals().second << "\n";

    vertex.get_r_sites().first  = r_site;
    vertex.get_r_sites().second = r_site;
  }

  /*
  template<class vertex_pair_type, class parameters_type, class concurrency_type, class H_interactions_type>//, class MOMS_type>
  static void set_vertex(vertex_pair_type& vertex, 
			 parameters_type&  parameters, 
			 concurrency_type& concurrency,
			 H_interactions_type& H_interactions,
			 //MOMS_type&        MOMS,
			 MC_integration_method_type MC_integration_method)
  {
    switch(MC_integration_method)
      {
      case CT_AUX:
	set_vertex_CT_AUX(vertex, parameters, concurrency, H_interactions);//MOMS);
	break;
      case PCM:
	set_vertex_PCM(vertex, parameters, concurrency, H_interactions);//, MOMS);
	break;
      default:
	throw std::logic_error(__FUNCTION__);
      }
  }

private:

  template<class vertex_pair_type, class parameters_type, class concurrency_type, class H_interactions_type>//, class MOMS_type>
  static void set_vertex_CT_AUX(vertex_pair_type& vertex, 
				parameters_type&  parameters, 
				concurrency_type& concurrency,
				H_interactions_type& H_interactions)
				//MOMS_type&        MOMS)
  {
    typedef typename parameters_type::electron_spin_domain_type electron_spin_domain_type;
    typedef typename parameters_type::DCA_cluster_type          DCA_cluster_type;
    typedef typename parameters_type::nu                        nu;
    typedef typename parameters_type::nu_nu_r_DCA               nu_nu_r_DCA;

    static int BANDS                                     = parameters.get_interacting_bands().size();
//     static FUNC_LIB::function<double, nu_nu_r_DCA>& H_interactions = MOMS.H_interactions;

    do
      {
	int band_ind_1 = concurrency.get_random_number()*BANDS;
	int band_ind_2 = concurrency.get_random_number()*BANDS;
      
	vertex.get_bands().first  = parameters.get_interacting_bands()[band_ind_1];
	vertex.get_bands().second = parameters.get_interacting_bands()[band_ind_2];
      
	vertex.get_e_spins().first = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	if(vertex.get_bands().first == vertex.get_bands().second)
	  {
	    if(vertex.get_e_spins().first == e_UP)
	      vertex.get_e_spins().second = e_DN;
	    else
	      vertex.get_e_spins().second = e_UP;
	  }
	else
	  vertex.get_e_spins().second = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	vertex.get_spin_orbitals().first  = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().first , vertex.get_e_spins().first); 
	vertex.get_spin_orbitals().second = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().second, vertex.get_e_spins().second); 
      }
    while(std::fabs(H_interactions(vertex.get_spin_orbitals().first, vertex.get_spin_orbitals().second, 0)) < 1.e-3 );

    int r_site = int(DCA_cluster_type::get_cluster_size()*concurrency.get_random_number());

    //   cout << r_site << "\t" << vertex.get_spin_orbitals().first << "\t" << vertex.get_spin_orbitals().second << "\n";

    vertex.get_r_sites().first  = r_site;
    vertex.get_r_sites().second = r_site;
  }

  template<class vertex_pair_type, class parameters_type, class concurrency_type, class H_interactions_type>
  static void set_vertex_PCM(vertex_pair_type&    vertex, 
			     parameters_type&     parameters, 
			     concurrency_type&    concurrency,
			     H_interactions_type& H_interactions)
  {
    typedef typename parameters_type::electron_spin_domain_type electron_spin_domain_type;
    typedef typename parameters_type::PCM_cluster_type          PCM_cluster_type;
    typedef typename parameters_type::nu                        nu;
    typedef typename parameters_type::nu_nu_r_PCM               nu_nu_r_PCM;

    int BANDS = parameters.get_interacting_bands().size();
//     static FUNC_LIB::function<double, nu_nu_r_PCM>& H_interactions = MOMS.H_interactions;

    do
      {
	int band_ind_1 = concurrency.get_random_number()*BANDS;
	int band_ind_2 = concurrency.get_random_number()*BANDS;
      
	vertex.get_bands().first  = parameters.get_interacting_bands()[band_ind_1];
	vertex.get_bands().second = parameters.get_interacting_bands()[band_ind_2];
      
	vertex.get_e_spins().first = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	if(vertex.get_bands().first == vertex.get_bands().second)
	  {
	    if(vertex.get_e_spins().first == e_UP)
	      vertex.get_e_spins().second = e_DN;
	    else
	      vertex.get_e_spins().second = e_UP;
	  }
	else
	  vertex.get_e_spins().second = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];
      
	vertex.get_spin_orbitals().first  = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().first , vertex.get_e_spins().first); 
	vertex.get_spin_orbitals().second = QMC::convert<int, nu>::spin_orbital(vertex.get_bands().second, vertex.get_e_spins().second); 
      
      }
    while(std::fabs(H_interactions(vertex.get_spin_orbitals().first, vertex.get_spin_orbitals().second, 0)) < 1.e-3 );

    int r_site = int(PCM_cluster_type::get_cluster_size()*concurrency.get_random_number());

    //    cout << r_site << "\t" << vertex.get_spin_orbitals().first << "\t" << vertex.get_spin_orbitals().second << "\n";

    vertex.get_r_sites().first  = r_site;
    vertex.get_r_sites().second = r_site;
  }
  */
};

#endif
