//-*-C++-*-

/*! 
 *  author Peter Staar
 *  
 */

#ifndef cached_auxilery_field_values_H
#define cached_auxilery_field_values_H

namespace QMC {
  
  class CV
  {
#include "type_definitions.h" 

  public:

    static int    spin_orbital(int band, e_spin_states_type e_spin);

    static double gamma(int spin_orbital_1, int spin_orbital_2, int site=0);    

    static double exp_V(int spin_orbital_1, int spin_orbital_2, 
			HS_spin_states_type HS_spin, HS_field_sign_type HS_field_sign, int site=0);

    static double exp_delta_V(int spin_orbital_1, 
			      int spin_orbital_2, 
			      HS_spin_states_type HS_spin_1,
			      HS_spin_states_type HS_spin_2,
			      HS_field_sign_type HS_field_sign, 
			      int site=0);

    template<class parameter_type, class MOMS_type>
    static void initialize(parameter_type& parameters, MOMS_type& MOMS);

  private:
    
    static double BETA;
    static double K_CT_AUX;
    static double BANDS;
    static double FULL_CLUSTER_SIZE;
    static double CORRELATED_ORBITALS;

    static function<double, nu_nu_r_DCA>                 H_interaction;

    static function<int,    nu>&                         intitialize_spin_orbital();
    static function<double, nu_nu_r_DCA>&                initialize_gamma(); 
    static function<double, nu_nu_HS_s_HS_f_r_DCA>&      initialize_exp_V(); 
    static function<double, nu_nu_HS_s_HS_s_HS_f_r_DCA>& initialize_exp_delta_V();
  };

  double CV::BETA                = 0;
  double CV::K_CT_AUX            = 1;
  double CV::BANDS               = 1;    
  double CV::FULL_CLUSTER_SIZE   = 1;
  double CV::CORRELATED_ORBITALS = 0;

  function<double, CV::nu_nu_r_DCA> CV::H_interaction;
 
  double CV::gamma(int spin_orbital_1, int spin_orbital_2, int site)
  {
    static function<double, nu_nu_r_DCA>& gamma_function = initialize_gamma();
    return gamma_function(spin_orbital_1, spin_orbital_2, site);
  }
 
  double CV::exp_V(int spin_orbital_1, 
		   int spin_orbital_2, 
		   HS_spin_states_type HS_spin, 
		   HS_field_sign_type HS_field_sign, 
		   int site)
  {
    static function<double, nu_nu_HS_s_HS_f_r_DCA>& exp_V_function = initialize_exp_V();
    int HS_spin_ind  = do_cast<int>::execute(HS_spin);
    int HS_field_ind = do_cast<int>::execute(HS_field_sign);
    return exp_V_function(spin_orbital_1, spin_orbital_2, HS_spin_ind, HS_field_ind, site);
  } 

  double CV::exp_delta_V(int spin_orbital_1, 
			 int spin_orbital_2, 
			 HS_spin_states_type HS_spin_1,
			 HS_spin_states_type HS_spin_2,
			 HS_field_sign_type HS_field_sign, 
			 int site)
  {
    static function<double, nu_nu_HS_s_HS_s_HS_f_r_DCA>& exp_V_function = initialize_exp_delta_V();

    int HS_spin_1_ind = do_cast<int>::execute(HS_spin_1);
    int HS_spin_2_ind = do_cast<int>::execute(HS_spin_2);
    int HS_field_ind  = do_cast<int>::execute(HS_field_sign);

    return exp_V_function(spin_orbital_1, spin_orbital_2, HS_spin_1_ind, HS_spin_2_ind, HS_field_ind, site);
  } 


  template<class parameter_type, class MOMS_type>
  void CV::initialize(parameter_type& parameters, MOMS_type& MOMS)
  {
    BETA               = time_domain_type::beta;
    K_CT_AUX           = parameters.get_K_CT_AUX();
    BANDS              = electron_band_domain_type::get_size();
    FULL_CLUSTER_SIZE  = DCA_cluster_type::get_cluster_size();

    H_interaction      = MOMS.H_interactions;

    CORRELATED_ORBITALS = 0;

    for(int nu_ind_i=0; nu_ind_i<2*BANDS; nu_ind_i++)
      for(int nu_ind_j=0; nu_ind_j<2*BANDS; nu_ind_j++)
	for(int r=0; r<FULL_CLUSTER_SIZE; r++)
	  fabs(H_interaction(nu_ind_i, nu_ind_j, r)) > 1.e-3 ? CORRELATED_ORBITALS++ : CORRELATED_ORBITALS ;

    CORRELATED_ORBITALS = CORRELATED_ORBITALS/2.; 
  }

  function<double, CV::nu_nu_r_DCA>&      CV::initialize_gamma()
  {
    //cout << __FUNCTION__ << endl;

    static function<double, nu_nu_r_DCA> gamma_function;
    
    for(int nu_ind_i=0; nu_ind_i<2*BANDS; nu_ind_i++)
      {
	for(int nu_ind_j=0; nu_ind_j<2*BANDS; nu_ind_j++)
	  {
	    for(int r=0; r<FULL_CLUSTER_SIZE; r++)
	      {
		double U_i_j_r = H_interaction(nu_ind_i, nu_ind_j, r);

		double coshgamma=1.+U_i_j_r*BETA*FULL_CLUSTER_SIZE*/*BANDS*/CORRELATED_ORBITALS/(2.*K_CT_AUX); 
		
		gamma_function(nu_ind_i, nu_ind_j, r) = acosh(coshgamma);

		//cout << gamma_function(nu_ind_i, nu_ind_j, r) << endl;
	      }
	  }
      }

    return gamma_function;
  }

  function<double, CV::nu_nu_HS_s_HS_f_r_DCA>& CV::initialize_exp_V()
  {
    //cout << __FUNCTION__ << endl;

    static function<double, nu_nu_HS_s_HS_f_r_DCA> exp_V_function;

    for(int nu_ind_i=0; nu_ind_i<2*BANDS; nu_ind_i++)
      {
	for(int nu_ind_j=0; nu_ind_j<2*BANDS; nu_ind_j++)
	  {
	    for(int HS_spin_ind=0; HS_spin_ind<3; HS_spin_ind++)
	      {
		for(int HS_field_ind=0; HS_field_ind<2; HS_field_ind++)
		  {
		    for(int r=0; r<FULL_CLUSTER_SIZE; r++)
		      {
			HS_spin_states_type HS_spin  = HS_spin_domain_type::get_elements()[HS_spin_ind];
			HS_field_sign_type  HS_field = HS_field_sign_domain_type::get_elements()[HS_field_ind];

			exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind, HS_field_ind, r) = 
			  std::exp(-gamma(nu_ind_i, nu_ind_j, r)*HS_spin*HS_field);

// 			cout << -gamma(nu_ind_i, nu_ind_j, r)*HS_spin*HS_field 
// 			     <<  "\t --> \t"
// 			     << exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind, HS_field_ind, r) 
// 			     << endl;
		      }
		  }
	      }
	  }
      }

    return exp_V_function;
  }

  function<double, CV::nu_nu_HS_s_HS_s_HS_f_r_DCA>& CV::initialize_exp_delta_V()
  {
    static function<double, nu_nu_HS_s_HS_s_HS_f_r_DCA> exp_delta_V_function;

    for(int nu_ind_i=0; nu_ind_i<2*BANDS; nu_ind_i++)
      {
	for(int nu_ind_j=0; nu_ind_j<2*BANDS; nu_ind_j++)
	  {
	    for(int HS_spin_1_ind=0; HS_spin_1_ind<3; HS_spin_1_ind++)
	      {
		for(int HS_spin_2_ind=0; HS_spin_2_ind<3; HS_spin_2_ind++)
		  {
		    for(int HS_field_ind=0; HS_field_ind<2; HS_field_ind++)
		      {
			for(int r=0; r<FULL_CLUSTER_SIZE; r++)
			  {
			    HS_spin_states_type HS_spin_1  = HS_spin_domain_type::get_elements()[HS_spin_1_ind];
			    HS_spin_states_type HS_spin_2  = HS_spin_domain_type::get_elements()[HS_spin_2_ind];

			    HS_field_sign_type  HS_field = HS_field_sign_domain_type::get_elements()[HS_field_ind];
			    
			    exp_delta_V_function(nu_ind_i, nu_ind_j, HS_spin_1_ind, HS_spin_2_ind, HS_field_ind, r) = 
			      std::exp(-gamma(nu_ind_i, nu_ind_j, r)*(HS_spin_1-HS_spin_2)*HS_field);
			  }
		      }
		  }
	      }
	  }
      }

    return exp_delta_V_function;
  }

}

#endif
