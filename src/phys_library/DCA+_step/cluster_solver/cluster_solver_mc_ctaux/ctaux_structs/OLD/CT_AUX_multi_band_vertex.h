//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file CT_AUX_multi_band_vertex.h  
 *
 */

#ifndef CT_AUX_multi_band_vertex_H
#define CT_AUX_multi_band_vertex_H


namespace dca {

  /*  template<class parameters_type, class concurrency_type>
  class vertex_couple
  {
#include "type_definitions.h" 

    typedef vertex_couple<parameters_type, concurrency_type> this_type;

  public:

    vertex_couple(parameters_type&   parameters_ref, concurrency_type&  concurrency_ref):
      parameters(parameters_ref),
      concurrency(concurrency_ref),
      interacting_bands(parameters_ref.get_interacting_bands()),
      BANDS(interacting_bands.size()),
      bands(0,0),
      e_spins(e_DN, e_UP),
      spin_orbitals(0,1),
      r_sites(0,0)
    {}

    ~vertex_couple()
    {}

    this_type& operator=(this_type& other_vertex_couple) // --> necessary for push_back
    {
      bands         = other_vertex_couple.get_bands();
      e_spins       = other_vertex_couple.get_e_spins();
      spin_orbitals = other_vertex_couple.get_spin_orbitals();
      r_sites       = other_vertex_couple.get_r_sites();

      return *this;
    }

    void set_random()
    {
      bands.first  = concurrency.get_random_number()*BANDS;
      bands.second = concurrency.get_random_number()*BANDS;
      
      e_spins.first = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];

      if(bands.first == bands.second)
	{
	  if(e_spins.first == e_UP)
	    e_spins.second = e_DN;
	  else
	    e_spins.second = e_UP;
	}
      else
	e_spins.second = electron_spin_domain_type::get_elements()[int(concurrency.get_random_number()*2.)];

      spin_orbitals.first  = CV::spin_orbital(bands.first , e_spins.first); 
      spin_orbitals.second = CV::spin_orbital(bands.second, e_spins.second); 

      r_sites.first  = int(DCA_cluster::get_cluster_size()*concurrency.get_random_number());
      r_sites.second = r_sites.first;

      assert(spin_orbitals.first != spin_orbitals.second);
    }

    
    std::pair<int,int>& get_bands()
    {
      return bands;
    }

    std::pair<e_spin_states_type, e_spin_states_type>& get_e_spins()
    {
      return e_spins;
    }

    std::pair<int,int>& get_spin_orbitals()
    {
      return spin_orbitals;
    }

    std::pair<int,int>& get_r_sites()
    {
      return r_sites;
    }

  private:

    parameters_type&    parameters;
    concurrency_type&   concurrency;

    std::vector<int>&  interacting_bands;
    int                BANDS;

    // pairs over the HS_field_signs
    std::pair<int,int>                                bands;
    std::pair<e_spin_states_type, e_spin_states_type> e_spins;

    std::pair<int,int>                                spin_orbitals;
    std::pair<int,int>                                r_sites;
  };
  */















  template<class parameters_type, class MOMS_type, class concurrency_type>
  class CT_AUX_multi_band_vertex
  {
#include "type_definitions.h" 

  public:
    
    typedef HS_spin_states_type                                                    spin_state_type;
    typedef CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type> this_type;
    typedef vertex_pair<parameters_type, concurrency_type>                       vertex_pair_type;

  public:

    CT_AUX_multi_band_vertex(parameters_type&   parameters_ref,
			      MOMS_type&         MOMS_ref,
			      concurrency_type&  concurrency_ref);
    
    ~CT_AUX_multi_band_vertex();

    this_type& operator=(this_type& vertex); // --> necessary for push_back

    void set_random();

    void set_random_interacting();
    void set_random_noninteracting();

    vertex_pair_type&  get_vertices();
    std::pair<int,int>&   get_bands();
    std::pair<e_spin_states_type, e_spin_states_type>&  get_e_spins();
    std::pair<int,int>&   get_sites();
    std::pair<int,int>&   get_spin_orbitals();
    

    HS_spin_states_type& get_HS_spin();
    int&                 get_delta_r();
    double&              get_tau();

  private:

    parameters_type&    parameters;
    MOMS_type&          MOMS;
    concurrency_type&   concurrency;

    vertex_pair_type  vertices;

    HS_spin_states_type HS_spin;
    int                 delta_r;
    double              tau;
    
  };

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::CT_AUX_multi_band_vertex(parameters_type&   parameters_ref,
												   MOMS_type&         MOMS_ref,
												   concurrency_type&  concurrency_ref):
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    concurrency(concurrency_ref),
    
    vertices(parameters_ref, concurrency_ref),
    HS_spin(HS_ZERO),
    delta_r(0),
    tau(0.)
  {}

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::~CT_AUX_multi_band_vertex()
  {}

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>& 
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::operator=(this_type& other_vertex) // --> necessary for push_back
  {
    vertices = other_vertex.get_vertices();
    HS_spin  = other_vertex.get_HS_spin();
    delta_r  = other_vertex.get_delta_r();
    tau      = other_vertex.get_tau();
    
    return *this;
  }

  
  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random()
  {
    vertices.set_random();

    double draw = concurrency.get_random_number();

    if(draw > 2./3.)
      HS_spin = HS_UP;
    else
      if(draw > 1./3.)
	HS_spin = HS_DN;
      else
	HS_spin = HS_ZERO;

    delta_r = DCA_r_cluster_type::subtract(vertices.get_r_sites().first, vertices.get_r_sites().second);

    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random_interacting()
  {
    vertices.set_random();

    double draw = concurrency.get_random_number();

    if(draw > 1./2.)
      HS_spin = HS_UP;
    else
      HS_spin = HS_DN;

    delta_r = DCA_r_cluster_type::subtract(vertices.get_r_sites().first, vertices.get_r_sites().second);

    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random_noninteracting()
  {
    vertices.set_random();

    HS_spin = HS_ZERO;

    delta_r = DCA_r_cluster_type::subtract(vertices.get_r_sites().first, vertices.get_r_sites().second);

    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  vertex_pair<parameters_type, concurrency_type>& 
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_vertices()
  {
    return vertices;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  std::pair<int,int>&   CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_bands()
  {
    return vertices.get_bands();
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  std::pair<e_spin_states_type, e_spin_states_type>&  
  CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_e_spins()
  {
    return vertices.get_e_spins();
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  std::pair<int,int>& CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_sites()
  {
    return vertices.get_r_sites();
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  std::pair<int,int>& CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_spin_orbitals()
  {
    return vertices.get_spin_orbitals();
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  HS_spin_states_type& CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_HS_spin()
  {
    return HS_spin;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  int& CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_delta_r()
  {
    return delta_r;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  double& CT_AUX_multi_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_tau()
  {
    return tau;
  }
  






} 

#endif

/*@}*/
