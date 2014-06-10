//-*-C++-*-

#ifndef SS_CT_HYBRIDIZATION_WALKER_H
#define SS_CT_HYBRIDIZATION_WALKER_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \defgroup CT-AUX-WALKER
     *  \ingroup  SS-CT-HYB
     */

    /*!
     *  \defgroup STRUCTURES
     *  \ingroup  SS-CT-HYB
     */

    /*!
     *  \ingroup SS-CT-HYB
     *
     *  \brief   This class organizes the MC-walker in the SS CT-HYB QMC
     *  \author  Bart Ydens, Peter Staar, Andrei Plamada
     *  \version 1.0
     */
    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    class MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>
    {
#include "type_definitions.h"

      typedef typename parameters_type::rng_type rng_type;

      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::profiler_type    profiler_type;
      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::concurrency_type concurrency_type;

    public:

      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::vertex_vertex_matrix_type vertex_vertex_matrix_type;
      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::configuration_type configuration_type;

      typedef function<vertex_vertex_matrix_type, nu> M_matrix_type;

      typedef ss_hybridization_solver_routines<parameters_type, MOMS_type>                               ss_hybridization_solver_routines_type;
      typedef ss_hybridization_walker_routines<parameters_type, MOMS_type, configuration_type, rng_type> ss_hybridization_walker_routines_type;

      typedef full_line_tools    <ss_hybridization_walker_routines_type>     full_line_tools_t;
      typedef anti_segment_tools <ss_hybridization_walker_routines_type>  anti_segment_tools_t;
      typedef segment_tools      <ss_hybridization_walker_routines_type>       segment_tools_t;
      typedef shift_segment_tools<ss_hybridization_walker_routines_type> shift_segment_tools_t;
      typedef swap_segment_tools <ss_hybridization_walker_routines_type>  swap_segment_tools_t;

    public:

      MC_walker(parameters_type& parameters_ref,
                MOMS_type&       MOMS_ref,
                rng_type&        rng_ref,
                int              id=0);

      ~MC_walker();

      /*!
       *  \brief Initializes the configuration and sets \f$\mu_i = \frac12 \sum_j \frac{U_{ij}+U_{ji}}{2}\f$.
       */
      void  initialize();//function<double, nu> mu_DC);

      /*!
       *  \brief Returns if the configuration has gone through a warm-up sweep.
       */
      bool& is_thermalized() { return thermalized; }

      /*!
       *  \brief Goes through an integration sweep. Do N insertion, removal, shift or swap steps.
       *  The insertion and removal steps include the insertion or removal of a full line, a segment or an anti-segment.
       */
      void  do_sweep();

      /*!
       *  \brief Goes through an integration step. Do N insertion, removal, shift or swap steps.
       *  The insertion and removal steps include the insertion or removal of a full line, a segment or an anti-segment.
       */
      void do_step();

      /*!
       *  \brief Returns the QMC sign.
       */
      double get_sign() { return sign; }

      /*!
       *  \brief Returns the current configuration.
       */
      configuration_type& get_configuration() { return configuration; }

      /*!
       *  \brief Returns the current inverse hybridization matrix \f$M = F_r(t)^{-1}\f$.
       */
      M_matrix_type& get_M_matrices() {return M; }

      /*!
       *  \brief Returns the hybridization_tools object
       */
      ss_hybridization_walker_routines_type& get_ss_hybridization_walker_routines() { return ss_hybridization_walker_routines_obj; }

      /*!
       *  \brief Print the hybridization functions \f$F_k(w)\f$, \f$F_k(t)\f$ and \f$F_r(t)\f$.
       */
      template<class stream_type>
      void to_JSON(stream_type& ss);

    private:

      void test_interpolation();

      int  get_random_interacting_flavor();

      //       void construct_F_k_w();
      //       void construct_F_r_t();
	void do_insert_remove             (int j);
      void insert_or_remove_full_line   (int j);
      void insert_or_remove_anti_segment(int j);
      void insert_or_remove_segment     (int j);
      void shift_segment                (int j);
      void swap_random_orbitals         ();

      //       // high-frquency FT ...
      //       void compute_moments (function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      //       void add_moments     (function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      //       void subtract_moments(function<std::complex<double>, nu_nu_k_DCA_w>& f_source);
      //       void compensate_for_moments(function<std::complex<double>, nu_nu_k_DCA_t>& f_source);

    private:

      parameters_type&  parameters;
      MOMS_type&        MOMS;
      concurrency_type& concurrency;

      rng_type&         rng;

      int               thread_id;

      configuration_type configuration;

      ss_hybridization_solver_routines_type ss_hybridization_solver_routines_obj;
      ss_hybridization_walker_routines_type ss_hybridization_walker_routines_obj;

      function<double, nu>&            mu;
      function<double, nu_nu_r_DCA_t>& F_r_t;

      full_line_tools_t     full_line_tools_obj;
      anti_segment_tools_t  anti_segment_tools_obj;
      segment_tools_t       segment_tools_obj;
      shift_segment_tools_t shift_segment_tools_obj;
      swap_segment_tools_t  swap_segment_tools_obj;

      //std::vector<bool> is_interacting_band;

      /*
        function<std::complex<double>, nu_nu_k_DCA_w> F_k_w;
        function<std::complex<double>, nu_nu_k_DCA_t> F_k_t;
        function<double              , nu_nu_r_DCA_t> F_r_t;

        function<double                   , nu> mu;
      */

      function<vertex_vertex_matrix_type, nu> M;

      //    public:
      /*
        function<double, nu> mu_HALF;

        function<double, nu> a0;
        function<double, nu> a1;
      */

      bool    thermalized;

      double  sign;

      size_t nb_updates;
      size_t nb_successfull_updates;

      double percentage_steps ;
      double percentage_shifts;
      double percentage_swaps ;
      
      double total;
      
      double p_0;
      double p_1;
    };

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::MC_walker(parameters_type& parameters_ref,
                                                                          MOMS_type&       MOMS_ref,
                                                                          rng_type&        rng_ref,
                                                                          int              id):
      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      rng(rng_ref),

      thread_id(id),

      configuration(),

      ss_hybridization_solver_routines_obj(parameters, MOMS),
      ss_hybridization_walker_routines_obj(parameters, MOMS, configuration, rng),

      mu   (ss_hybridization_solver_routines_obj.get_mu()),
      F_r_t(ss_hybridization_solver_routines_obj.get_F_r_t()),

      full_line_tools_obj    (ss_hybridization_walker_routines_obj),
      anti_segment_tools_obj (ss_hybridization_walker_routines_obj),
      segment_tools_obj      (ss_hybridization_walker_routines_obj),
      shift_segment_tools_obj(ss_hybridization_walker_routines_obj),
      swap_segment_tools_obj (ss_hybridization_walker_routines_obj),

      /*
        F_k_w("cluster_hybridization_F_k_w"),
        F_k_t("cluster_hybridization_F_k_t"),
        F_r_t("cluster_hybridization_F_r_t"),

        mu("orbital_dependent_chemical_potential"),
      */

      M("M-matrices"),

      thermalized(false),

      sign(1)
    {
      //FLAVORS = b::dmn_size()*s::dmn_size();
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::~MC_walker()
    {
      if(concurrency.id()==0 and thread_id==0)
	{
	  stringstream ss;
	  ss << "\n\n\t\t walker died --> nb_successfull_updates/nb_updates : " << double(nb_successfull_updates)/double(nb_updates) << "\n\n";
	  cout << ss.str();
	}
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::initialize()
    {
      ss_hybridization_solver_routines_obj.initialize_functions();

      ss_hybridization_walker_routines_obj.initialize_akima_coefficients(F_r_t);

      //test_interpolation();

      {
	sign = 1;

	nb_updates             = 0;
	nb_successfull_updates = 0;
      }

      {
	percentage_steps  = parameters.get_steps_per_sweep();
	percentage_shifts = parameters.get_shifts_per_sweep();
	
	total = (percentage_steps+percentage_shifts);
	
	p_0 = (percentage_steps + 0                )/total;
	p_1 = (percentage_steps + percentage_shifts)/total;
      }

      {
	configuration.initialize();
	
	is_thermalized() = false;
	
	for(int i=0; i<M.size(); i++){
	  M(i).get_current_size() = std::pair<int,int>(0,0);
	  
	  //M(i).print_fingerprint();
	}
      }
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::test_interpolation()
    {
      cout << __FUNCTION__ << endl;

      SHOW::execute_on_bands(F_r_t);

      {
        double beta = parameters.get_beta();

        int coor[2];

        for(int b_i=0; b_i<b::dmn_size(); b_i++){
          for(int s_i=0; s_i<s::dmn_size(); s_i++){

            coor[0] = b_i;
            coor[1] = s_i;

            std::vector<double> x(0);
            std::vector<double> y(0);

            for(int t_i=0; t_i<10*t::dmn_size(); t_i++){

              double t_val = -beta + 2*beta/(10*t::dmn_size()-1)*t_i;

              if(t_i==0)
                t_val += 1.e-6;

              if(t_i==10*t::dmn_size()-1)
                t_val -= 1.e-6;

              double F_val = ss_hybridization_walker_routines_obj.interpolate_F(coor, t_val, F_r_t);

              x.push_back(t_val);
              y.push_back(F_val);
            }

            SHOW::plot_points(x, y);
          }
        }
      }

      throw std::logic_error(__FUNCTION__);
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_sweep()
    {
      double factor = 1.;
      if(thermalized)
        factor = parameters.get_number_of_sweeps_per_measurement();

      int nr_of_segments = std::max(16, configuration.size());

      double ratio  = double(nb_successfull_updates+1)/double(nb_updates+1);
      int    factor2 = std::max(1, int(std::ceil(1.0/ratio)));

      int nb_steps = nr_of_segments*factor*factor2;

      for(int l=0; l<nb_steps; l++)
        {
          if(true)
            {
              do_step();
            }
	  /*
          else
            {
              {
                int so_ind = get_random_interacting_flavor();

                if(configuration.get_vertices(so_ind).size() == 0)
                  insert_or_remove_full_line(so_ind);
              }

              {
                {
                  int so_ind = get_random_interacting_flavor();

                  insert_or_remove_anti_segment(so_ind);
                }

                {
                  int so_ind = get_random_interacting_flavor();

                  if(!configuration.get_full_line(so_ind))
                    insert_or_remove_segment(so_ind);
                }
              }
            }
	  */
        }
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    int MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::get_random_interacting_flavor()
    {
      int spin     = s::dmn_size()*rng.get_random_number();
      int int_band = parameters.get_interacting_bands().size()*rng.get_random_number();
      
      return parameters.get_interacting_bands()[int_band] + spin*b::dmn_size();
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_step()
    {
      if(true)
        {
	  /*
          double percentage_steps  = parameters.get_steps_per_sweep();
          double percentage_shifts = parameters.get_shifts_per_sweep();
          double percentage_swaps  = 0;//parameters.get_swaps_per_sweep();

          double total = (percentage_steps+percentage_shifts+percentage_swaps);

          double p_0 = (percentage_steps + 0                 + 0.0             )/total;
          double p_1 = (percentage_steps + percentage_shifts + 0.0             )/total;
	  */




          double p = rng.get_random_number();
	  
	  int so_ind = get_random_interacting_flavor();

          if( p < p_0)
            {
		do_insert_remove(so_ind);
            }

          if( p_0 < p and p < p_1)
            {
              shift_segment(so_ind);
            }

//           if( p_1 < p)
//             {
//               swap_random_orbitals();
//             }
        }
      /*
      else
        {
          {
            int so_ind = get_random_interacting_flavor();

            if(configuration.get_vertices(so_ind).size() == 0)
              insert_or_remove_full_line(so_ind);
          }

          {
            {
              int so_ind = get_random_interacting_flavor();

              insert_or_remove_anti_segment(so_ind);
            }

            {
              int so_ind = get_random_interacting_flavor();

              if(!configuration.get_full_line(so_ind))
                insert_or_remove_segment(so_ind);
            }
          }

        }
      */

      /*
        {
        int j = get_random_interacting_flavor();

        if(configuration.get_vertices(j).size() == 0)
        insert_or_remove_full_line(j);
        }

        for(int k=0; k<parameters.get_steps_per_sweep(); k++)
        {
        int j = get_random_interacting_flavor();

        insert_or_remove_anti_segment(j);

        if(!configuration.get_full_line(j))
        insert_or_remove_segment(j);
        }
      */

//       if(concurrency.id()==0 and id==0)
// 	{
// 	  stringstream ss;
// 	  ss << "\t" << nb_successfull_updates << ", " << nb_updates << ", " << nb_successfull_updates/nb_updates << "\n\n";
// 	  cout << ss.str();
// 	}
    }


	template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
	void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::do_insert_remove(int so_ind)
	{
	    double rn=rng.get_random_number();

	    if(configuration.get_vertices(so_ind).size() == 0)
	    {
		bool succes;
		nb_updates += 1;

	    	if(configuration.get_full_line(so_ind))
		{
	    	    if(rn < 0.75)
	    		succes = full_line_tools_obj.insert_or_remove(so_ind, mu(so_ind));
	    	    else
			succes = anti_segment_tools_obj.insert_anti_segment(so_ind, mu(so_ind), sign, M, F_r_t);
		}
		else
		{
		    if(rn<0.75)
	    		succes = full_line_tools_obj.insert_or_remove(so_ind, mu(so_ind));
		    else
			succes = segment_tools_obj.insert_segment(so_ind, mu(so_ind), sign, M, F_r_t);
		}
		nb_successfull_updates += succes? 1 : 0;
	    }
	    else
	    {
		if(rn < 0.5)
		{
		    insert_or_remove_anti_segment(so_ind);
		}
		else
		{
		    insert_or_remove_segment(so_ind);
		}
	    }
	}

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_full_line(int j)
    {
      nb_updates += 1;

      bool succes = full_line_tools_obj.insert_or_remove(j, mu(j));

      nb_successfull_updates += succes? 1 : 0;
      
      if(QMC_INTEGRATOR_BIT)
        ss_ct_hybridization_walker_bit::FULL_CHECK(configuration, M, F_r_t, ss_hybridization_walker_routines_obj);
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_anti_segment(int j)
    {
      nb_updates += 1;

      bool succes;
      if(rng.get_random_number() < 0.50000000)
        succes = anti_segment_tools_obj.insert_anti_segment(j, mu(j), sign, M, F_r_t);
      else
        succes = anti_segment_tools_obj.remove_anti_segment(j, mu(j), sign, M, F_r_t);

      nb_successfull_updates += succes? 1 : 0;

      if(QMC_INTEGRATOR_BIT)
        ss_ct_hybridization_walker_bit::FULL_CHECK(configuration, M, F_r_t, ss_hybridization_walker_routines_obj);
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::insert_or_remove_segment(int j)
    {
      nb_updates += 1;

      bool succes;
      if(rng.get_random_number() < 0.50000000)
        succes = segment_tools_obj.insert_segment(j, mu(j), sign, M, F_r_t);
      else
        succes = segment_tools_obj.remove_segment(j, mu(j), sign, M, F_r_t);

      nb_successfull_updates += succes? 1 : 0;

      if(QMC_INTEGRATOR_BIT)
        ss_ct_hybridization_walker_bit::FULL_CHECK(configuration, M, F_r_t, ss_hybridization_walker_routines_obj);
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::shift_segment(int j)
    {
      nb_updates += 1;

      bool succes;
      if(rng.get_random_number() < 0.50000000)
	  succes= shift_segment_tools_obj.shift_segment_start_vertex(j, mu(j), sign, M, F_r_t);
      else
	  succes= shift_segment_tools_obj.shift_segment_end_vertex(j, mu(j), sign, M, F_r_t);

      nb_successfull_updates += succes? 1 : 0;

      if(QMC_INTEGRATOR_BIT)
        ss_ct_hybridization_walker_bit::FULL_CHECK(configuration, M, F_r_t, ss_hybridization_walker_routines_obj);
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::swap_random_orbitals()
    {
      nb_updates += 1;

      int i = get_random_interacting_flavor();
      int j = get_random_interacting_flavor();

      bool succes= swap_segment_tools_obj.swap_orbitals(i, j, mu, sign, M, F_r_t);
      nb_successfull_updates += succes? 1 : 0;

      if(QMC_INTEGRATOR_BIT)
        ss_ct_hybridization_walker_bit::FULL_CHECK(configuration, M, F_r_t, ss_hybridization_walker_routines_obj);
    }

    /*!
     *
     *   \brief We follow the convention that F(\omega) = omega*i+mu - G^{-1}_xc(\omega) (Notice that F behaves as H_0).
     */
    /*
      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::construct_F_k_w()
      {
      if(thread_id==0){
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
      std::complex<double> S = MOMS.Sigma(nu_ind, nu_ind, 0, w_ind);

      F_k_w(nu_ind, nu_ind, 0, w_ind) = z-(1./G+S);
      }
      }

      if(thread_id==0)
      SHOW::execute_on_bands(F_k_w);
      }
    */

    /*
      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::construct_F_r_t()
      {
      //cout << "\n\t\t implement " << __FUNCTION__ << " properly !!\n\n";
      //throw std::logic_error(__FUNCTION__);

      compute_moments(F_k_w);

      subtract_moments(F_k_w);

      //FT<w,t>::execute(F_k_w, F_k_t);
      MATH_ALGORITHMS::TRANSFORM<w, t>::execute(F_k_w, F_k_t);

      add_moments(F_k_w);

      compensate_for_moments(F_k_t);

      //FT<k_DCA,r_DCA>::execute(F_k_t,F_r_t);
      MATH_ALGORITHMS::TRANSFORM<k_DCA, r_DCA>::execute(F_k_t,F_r_t);

      //       for(int i=0; i<F_r_t.size(); i++)
      //         F_r_t(i) = abs(F_r_t(i));
      }

      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compute_moments(function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
      {
      int    w_ind = 0;
      double w_val = w::get_elements()[w_ind];

      for(int nu_ind=0;nu_ind<nu::dmn_size();nu_ind++)
      {
      a0(nu_ind) = real( F_k_w(nu_ind,nu_ind,0,w_ind));
      a1(nu_ind) = imag(-F_k_w(nu_ind,nu_ind,0,w_ind)*w_val);
      }
      }

      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::subtract_moments(function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
      {

      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){

      std::complex<double> w_val = std::complex<double>(0, w::get_elements()[w_ind]);

      for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++)
      f_source(nu_ind, nu_ind, 0, w_ind) -= (a0(nu_ind) + a1(nu_ind)/w_val);
      }

      if(thread_id==0)
      SHOW::execute(f_source);
      }

      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::add_moments(function<std::complex<double>, nu_nu_k_DCA_w>& f_source)
      {
      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){

      std::complex<double> w_val = std::complex<double>(0, w::get_elements()[w_ind]);

      for(int nu_ind=0; nu_ind<nu::dmn_size(); nu_ind++)
      f_source(nu_ind, nu_ind, 0, w_ind) += (a0(nu_ind) + a1(nu_ind)/w_val);
      }

      if(thread_id==0)
      SHOW::execute(f_source);
      }

      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type>::compensate_for_moments(function<std::complex<double>, nu_nu_k_DCA_t>& f_source)
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
    */
  }

}

#endif

