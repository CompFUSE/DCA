//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file CT_AUX_single_band_vertex.h  
 *
 */

#ifndef CT_AUX_single_band_vertex_H
#define CT_AUX_single_band_vertex_H


namespace dca {

  template<class parameters_type, class MOMS_type, class concurrency_type>
  class CT_AUX_single_band_vertex
  {
#include "type_definitions.h" 

  public:
    
    typedef CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type> this_type;
    typedef HS_spin_states_type spin_state_type;

  public:

    CT_AUX_single_band_vertex(parameters_type&   parameters_ref,
			      MOMS_type&         MOMS_ref,
			      concurrency_type&  concurrency_ref);
    
    CT_AUX_single_band_vertex(parameters_type&   parameters_ref,
			      MOMS_type&         MOMS_ref,
			      concurrency_type&  concurrency_ref,
			      HS_spin_states_type spin_val,
			      int                site_val,
			      double             tau_val);

    ~CT_AUX_single_band_vertex();

    this_type& operator=(this_type& vertex); // --> necessary for push_back
    

    void set_random();

    void set_random_interacting();
    void set_random_noninteracting();

    HS_spin_states_type& get_spin();
    int&                 get_site();
    int&                 get_band();
    double&              get_tau();

  private:

    parameters_type&                   parameters;
    MOMS_type&                         MOMS;
    concurrency_type&                  concurrency;

    HS_spin_states_type spin;
    int                 site;
    static int          band;
    double              tau;
    
  };

  template<class parameters_type, class MOMS_type, class concurrency_type>
  int CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::band = 0;


  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::CT_AUX_single_band_vertex(parameters_type&   parameters_ref,
												     MOMS_type&         MOMS_ref,
												     concurrency_type&  concurrency_ref):
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    concurrency(concurrency_ref),
    spin(HS_ZERO),
    site(0),
    tau(0.)
  {  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::CT_AUX_single_band_vertex(parameters_type&   parameters_ref,
												     MOMS_type&         MOMS_ref,
												     concurrency_type&  concurrency_ref,
												     HS_spin_states_type spin_val,
												     int                site_val,
												     double             tau_val):
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    concurrency(concurrency_ref),
    spin(spin_val),
    site(site_val),
    tau(tau_val)
  {
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::~CT_AUX_single_band_vertex()
  {}

  template<class parameters_type, class MOMS_type, class concurrency_type>
  CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>& 
  CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::operator=(this_type& vertex) // --> necessary for push_back
  {
    spin = vertex.get_spin();
    site = vertex.get_site();
    tau = vertex.get_tau();
    band = vertex.get_band();
    
    return *this;
  }


  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random()
  {
    double draw = concurrency.get_random_number();

    if(draw > 2./3.)
      spin = HS_UP;
    else
      if(draw > 1./3.)
	spin = HS_DN;
      else
	spin = HS_ZERO;

    site = concurrency.get_random_number()*DCA_cluster::get_cluster_size();

    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random_interacting()
  {
    double draw = concurrency.get_random_number();

    if(draw > 1./2.)
      spin = HS_UP;
    else
      spin = HS_DN;

    site = concurrency.get_random_number()*DCA_cluster::get_cluster_size();
    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  void CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::set_random_noninteracting()
  {
    //cout << __FUNCTION__ << endl;

    spin = HS_ZERO;

    site = concurrency.get_random_number()*DCA_cluster::get_cluster_size();
    tau = concurrency.get_random_number()*time_domain_type::beta;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  HS_spin_states_type& CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_spin()
  {
    return spin;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  int& CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_site()
  {
    return site;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  int& CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_band()
  {
    return band;
  }

  template<class parameters_type, class MOMS_type, class concurrency_type>
  double& CT_AUX_single_band_vertex<parameters_type, MOMS_type, concurrency_type>::get_tau()
  {
    return tau;
  }







} 

#endif

/*@}*/
