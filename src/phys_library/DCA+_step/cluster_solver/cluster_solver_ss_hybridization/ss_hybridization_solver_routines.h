//-*-C++-*-

#ifndef SS_HYBRIDIZATION_SOLVER_ROUTINES_H
#define SS_HYBRIDIZATION_SOLVER_ROUTINES_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \class   hybridization_routines
     *  \ingroup SS_HYBRIDIZATION
     *
     *  \author  Bart Ydens, Peter Staar
     *  \brief   This class implements the help functions for the insertion and removal of (anti-)segments. The help functions
     *  include the calculation of the deteminant ratio and the computation of the new hybridization matrix, using sherman-morrison equations.
     *
     */
    template<typename parameters_t, typename MOMS_t>
    class ss_hybridization_solver_routines
    {
#include "type_definitions.h"

      const static bool SHOW_FUNCTIONS = false;

    public:

      typedef parameters_t    parameters_type;
      typedef MOMS_t          MOMS_type;

      typedef typename parameters_type::profiler_type    profiler_type;
      typedef typename parameters_type::concurrency_type concurrency_type;

    public:

      ss_hybridization_solver_routines(parameters_t& parameters_ref,
                                       MOMS_t&       MOMS_ref);

      ~ss_hybridization_solver_routines();

      void initialize();

      void initialize_interacting_band_vector();

      void initialize_functions();

      bool is_interacting_band(int b_ind);

      FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& get_F_k_w() { return F_k_w; }
      FUNC_LIB::function<             double , nu_nu_r_DCA_t>& get_F_r_t() { return F_r_t; }

      FUNC_LIB::function<double, nu>& get_mu()      { return mu;}
      FUNC_LIB::function<double, nu>& get_mu_HALF() { return mu_HALF;}

      FUNC_LIB::function<double, nu>& get_a0() { return a0; }
      FUNC_LIB::function<double, nu>& get_a1() { return a1; }

    private:

      void initialize_hybridization_function();
      void initialize_hybridization_chem_pot();

      void construct_F_k_w();
      void construct_F_r_t();

      // high-frquency FT ...
      void compute_moments       (FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      void add_moments           (FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      void subtract_moments      (FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      void compensate_for_moments(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_t>& f_source);

    private:

      parameters_type&  parameters;
      MOMS_type&        MOMS;
      concurrency_type& concurrency;

      std::vector<bool> is_interacting_band_vector;

      FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w> F_k_w;
      FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_t> F_k_t;
      FUNC_LIB::function<double              , nu_nu_r_DCA_t> F_r_t;

      FUNC_LIB::function<double, nu> mu;
      FUNC_LIB::function<double, nu> mu_HALF;

      FUNC_LIB::function<double, nu> a0;
      FUNC_LIB::function<double, nu> a1;
    };

    template<typename parameters_t, typename MOMS_t>
    ss_hybridization_solver_routines<parameters_t, MOMS_t>::ss_hybridization_solver_routines(parameters_t&    parameters_ref,
                                                                                             MOMS_t&          MOMS_ref):
      parameters(parameters_ref),
      MOMS      (MOMS_ref),
      concurrency(parameters.get_concurrency()),

      is_interacting_band_vector(b::dmn_size(), false),

      F_k_w("cluster_hybridization_F_k_w"),
      F_k_t("cluster_hybridization_F_k_t"),
      F_r_t("cluster_hybridization_F_r_t"),

      mu     ("orbital_dependent_chemical_potential"),
      mu_HALF("mu-half"),

      a0("a0"),
      a1("a1")
    {
      initialize();
    }

    template<typename parameters_t, typename MOMS_t>
    ss_hybridization_solver_routines<parameters_t, MOMS_t>::~ss_hybridization_solver_routines()
    {}

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::initialize()
    {
      if(SHOW_FUNCTIONS)
        cout << "\n\t " << __FUNCTION__ << " \n";

      initialize_interacting_band_vector();
    }

    template<typename parameters_t, typename MOMS_t>
    bool ss_hybridization_solver_routines<parameters_t, MOMS_t>::is_interacting_band(int b_ind)
    {
      assert(b_ind>-1 and b_ind<b::dmn_size());
      return is_interacting_band_vector[b_ind];
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::initialize_interacting_band_vector()
    {
      if(SHOW_FUNCTIONS)
        cout << "\n\t " << __FUNCTION__ << " \n";

      for(int b_i=0; b_i<b::dmn_size(); b_i++){
        is_interacting_band_vector[b_i] = false;

        for(int b_j=0; b_j<parameters.get_interacting_bands().size(); b_j++)
          if(parameters.get_interacting_bands()[b_j]==b_i)
            is_interacting_band_vector[b_i] = true;
      }
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::initialize_functions()
    {
      initialize_hybridization_function();
      initialize_hybridization_chem_pot();
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::initialize_hybridization_chem_pot()
    {
      for(int b_i=0; b_i<b::dmn_size(); b_i++){

        if(is_interacting_band(b_i)){

          for(int s_i=0; s_i<s::dmn_size(); s_i++){

            mu     (b_i,s_i) = 0;
            mu_HALF(b_i,s_i) = 0;

            for(int s_j=0; s_j<s::dmn_size(); s_j++)
              for(int b_j=0; b_j<b::dmn_size(); b_j++)
                if(is_interacting_band(b_j))
                  mu_HALF(b_i,s_i) += (1./2.)*(MOMS.H_interactions(b_j,s_j,b_i,s_i,0)+MOMS.H_interactions(b_i,s_i,b_j,s_j,0))*1./2.;

//             if(parameters.get_double_counting_method() == "constant-correction")
	    {
	      mu(b_i,s_i) += (parameters.get_chemical_potential()
			      // + parameters.get_double_counting_correction()
			      + mu_HALF(b_i,s_i)
			      - a0     (b_i,s_i));
	    }
//             else
//               {
//                 mu(b_i,s_i) += (parameters.get_chemical_potential()
//                                 + parameters.get_double_counting_correction()
//                                 + mu_HALF(b_i,s_i)
//                                 - a0     (b_i,s_i));
//               }
          }
        }
      }
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::initialize_hybridization_function()
    {
      if(SHOW_FUNCTIONS)
        cout << "\n\t " << __FUNCTION__ << " \n";

      construct_F_k_w();

      if(SHOW_FUNCTIONS)
        SHOW::execute_on_bands(F_k_w);

      if(SHOW_FUNCTIONS)
        cout << "\n\t construct_F_r_t \n";

      construct_F_r_t();

      if(SHOW_FUNCTIONS)
        SHOW::execute_on_bands(F_r_t);

      //assert(false);
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::construct_F_k_w()
    {
      if(SHOW_FUNCTIONS)
        cout << "\n\t " << __FUNCTION__ << " \n";

      if(SHOW_FUNCTIONS)
        {
          SHOW::execute_on_bands(MOMS.G_k_w);
          SHOW::execute_on_bands(MOMS.Sigma);
          SHOW::execute_on_bands(MOMS.G0_k_w_cluster_excluded);
        }

      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){

        double mu = parameters.get_chemical_potential();
        double wm = w::parameter_type::get_elements()[w_ind];

        std::complex<double> z(mu, wm);

        for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++){

          std::complex<double> G = MOMS.G_k_w(nu_ind, nu_ind, 0, w_ind);
          //std::complex<double> S = MOMS.Sigma(nu_ind, nu_ind, 0, w_ind);
	  std::complex<double> S = MOMS.Sigma_cluster(nu_ind, nu_ind, 0, w_ind);

	  std::complex<double> G0_xc_inv = 0;

	  if(false)
	    G0_xc_inv = 1./MOMS.G0_k_w_cluster_excluded(nu_ind,nu_ind, 0, w_ind);
	  else
	    G0_xc_inv = (1./G+S);

          //F_k_w(nu_ind, nu_ind, 0, w_ind) = z-(1./G+S);
	  F_k_w(nu_ind, nu_ind, 0, w_ind) = z-G0_xc_inv;
        }
      }

      if(SHOW_FUNCTIONS)
        SHOW::execute_on_bands(F_k_w);
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::construct_F_r_t()
    {
      compute_moments(F_k_w);

      subtract_moments(F_k_w);

      MATH_ALGORITHMS::TRANSFORM<w, t>::execute(F_k_w, F_k_t);

      add_moments(F_k_w);

      compensate_for_moments(F_k_t);

      MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(F_k_t, F_r_t);
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::compute_moments(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
    {
      int    w_ind = 0;
      double w_val = w::get_elements()[w_ind];

      for(int nu_ind=0;nu_ind<nu::dmn_size();nu_ind++)
        {
          a0(nu_ind) = real( F_k_w(nu_ind,nu_ind,0,w_ind));
          a1(nu_ind) = imag(-F_k_w(nu_ind,nu_ind,0,w_ind)*w_val);
        }
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::subtract_moments(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
    {

      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){

        std::complex<double> w_val = std::complex<double>(0, w::get_elements()[w_ind]);

        for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++)
          f_source(nu_ind, nu_ind, 0, w_ind) -= (a0(nu_ind) + a1(nu_ind)/w_val);
      }

      //       if(SHOW_FUNCTIONS)
      //         SHOW::execute_on_bands(f_source);
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::add_moments(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
    {
      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){

        std::complex<double> w_val = std::complex<double>(0, w::get_elements()[w_ind]);

        for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++)
          f_source(nu_ind, nu_ind, 0, w_ind) += (a0(nu_ind) + a1(nu_ind)/w_val);
      }

      //       if(SHOW_FUNCTIONS)
      //         SHOW::execute_on_bands(f_source);
    }

    template<typename parameters_t, typename MOMS_t>
    void ss_hybridization_solver_routines<parameters_t, MOMS_t>::compensate_for_moments(FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_t>& f_source)
    {
      for(int t_ind=0; t_ind<t::dmn_size(); t_ind++){

        double t_val = t::get_elements()[t_ind];

        for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++){
          if(t_val<0)
            f_source(nu_ind, nu_ind, 0, t_ind) += a1(nu_ind)/2.;
          else
            f_source(nu_ind, nu_ind, 0, t_ind) -= a1(nu_ind)/2.;
        }
      }
    }

  }
}

#endif
