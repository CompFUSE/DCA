//-*-C++-*-

#ifndef DCA_QMCI_CTAUX_SP_ACCUMULATOR_NFFT_H
#define DCA_QMCI_CTAUX_SP_ACCUMULATOR_NFFT_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \ingroup CT-AUX
     *
     *  \brief   This class measures the single-particle functions with an NFFT scheme in the CT-AUX QMC engine.
     *  \author  Peter Staar
     *  \author  Raffaele Solca'
     *  \version 1.0
     */
    template<class parameters_type, class MOMS_type>
    class MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>
    {
#include "type_definitions.h"

      typedef vertex_singleton vertex_singleton_type;

      typedef r_DCA r_dmn_t;
      typedef k_DCA k_dmn_t;

      typedef w                      w_dmn_t;
      typedef dmn_3<nu, nu, r_dmn_t> p_dmn_t;

      typedef b b_dmn_t;
      typedef s s_dmn_t;

      typedef typename parameters_type::profiler_type    profiler_type;
      typedef typename parameters_type::concurrency_type concurrency_type;

      //typedef typename parameters_type::MC_measurement_scalar_type scalar_type;
      typedef double scalar_type;

    public:

      MC_single_particle_accumulator(parameters_type& parameters_ref);

      ~MC_single_particle_accumulator();

      void initialize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w);

      void initialize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w,
		      function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w_squared);

      void finalize();

      void finalize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w,
		    function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w_squared);

      template<class configuration_type, class vertex_vertex_matrix_type>
      void accumulate_M_r_w(configuration_type&        configuration_e_spin,
                            vertex_vertex_matrix_type& M,
                            double                     sign,
                            e_spin_states              e_spin);

      template<class configuration_type>
      void accumulate_K_r_t(configuration_type&                           configuration,
                            function<double, dmn_4<nu, nu, r_dmn_t, t> >& K_r_t,
                            double                                        sign);


      int find_first_non_interacting_spin(std::vector<vertex_singleton_type>& configuration_e_spin);

      void compute_M_r_w       (function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w);
      void compute_M_r_w_square(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w_square);

    private:

      parameters_type&  parameters;
      concurrency_type& concurrency;

      MATH_ALGORITHMS::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_M_r_w_obj;
      MATH_ALGORITHMS::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_M_r_w_squared_obj;
    };

    template<class parameters_type, class MOMS_type>
    MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::MC_single_particle_accumulator(parameters_type& parameters_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      cached_nfft_1D_M_r_w_obj(),
      cached_nfft_1D_M_r_w_squared_obj()
    {}

    template<class parameters_type, class MOMS_type>
    MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::~MC_single_particle_accumulator()
    {}

    template<class parameters_type, class MOMS_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::initialize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w)
    {
      cached_nfft_1D_M_r_w_obj.initialize();

      for(int i=0; i<M_r_w.size(); i++)
        M_r_w(i) = 0;
    }

    template<class parameters_type, class MOMS_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::initialize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w,
                                                                                                     function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w_squared)
    {
      cached_nfft_1D_M_r_w_obj        .initialize();
      cached_nfft_1D_M_r_w_squared_obj.initialize();

      for(int i=0; i<M_r_w.size(); i++)
        M_r_w(i) = 0;

      for(int i=0; i<M_r_w_squared.size(); i++)
        M_r_w_squared(i) = 0;
    }

    template<class parameters_type, class MOMS_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::finalize()
    {}

    template<class parameters_type, class MOMS_type>
    template<class configuration_type, class vertex_vertex_matrix_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::accumulate_M_r_w(configuration_type&        configuration_e_spin,
                                                                                                           vertex_vertex_matrix_type& M,
                                                                                                           double                     sign,
                                                                                                           e_spin_states              e_spin)
    {
      int coor_nfft[5];

      scalar_type one_div_two_beta = 1./(2.*parameters.get_beta());

      //int spin_index = do_cast<int>::execute(e_spin);
      int spin_index = electron_spin_domain::to_coordinate(e_spin);

      int         r_ind, b_i, b_j, r_i, r_j;//, s_i, s_j;
      scalar_type t_i, t_j, delta_tau, scaled_tau, f_tau;

      int configuration_size = find_first_non_interacting_spin(configuration_e_spin);//configuration_e_spin.size();

      for(int j=0; j<configuration_size; j++)
        {
          vertex_singleton_type& configuration_e_spin_j = configuration_e_spin[j];

          b_j = configuration_e_spin_j.get_band();
          r_j = configuration_e_spin_j.get_r_site();
          t_j = configuration_e_spin_j.get_tau();

          for(int i=0; i<configuration_size; i++)
            {
              //if(configuration_e_spin[i].get_HS_spin() != HS_ZERO && configuration_e_spin[j].get_HS_spin() != HS_ZERO)
              vertex_singleton_type& configuration_e_spin_i = configuration_e_spin[i];

              b_i = configuration_e_spin_i.get_band();
              r_i = configuration_e_spin_i.get_r_site();
              t_i = configuration_e_spin_i.get_tau();

              //            r_ind = r_cluster_type::subtract(r_j, r_i);
              r_ind = r_DCA::parameter_type::subtract(r_j, r_i);

              delta_tau = t_i-t_j;

              coor_nfft[0] = b_i;
              coor_nfft[1] = spin_index;
              coor_nfft[2] = b_j;
              coor_nfft[3] = spin_index;
              coor_nfft[4] = r_ind;

              scaled_tau = delta_tau*one_div_two_beta;

              scaled_tau += i==j ? 1.e-16 : 0.;

              f_tau = M(i,j)*sign;

              assert(configuration_e_spin[i].get_HS_spin() != HS_ZERO && configuration_e_spin[j].get_HS_spin() != HS_ZERO);

              cached_nfft_1D_M_r_w_obj        .accumulate_at(coor_nfft, scaled_tau, f_tau);
              cached_nfft_1D_M_r_w_squared_obj.accumulate_at(coor_nfft, scaled_tau, f_tau*f_tau);
            }
        }
    }

    template<class parameters_type, class MOMS_type>
    template<class configuration_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::accumulate_K_r_t(configuration_type&                           configuration,
                                                                                                           function<double, dmn_4<nu, nu, r_dmn_t, t> >& K_r_t,
                                                                                                           double                                        sign)
    {
      // for next generation solver !!

      /*
        typedef vertex_singleton<base_cluster_type> vertex_singleton_type;

        double delta_tau = 2.*parameters.get_beta()/(t::dmn_size()-2);

        for(int j=0; j<configuration.size(); j++)
        {
        vertex_singleton_type& vertex = configuration[j];

        //      int b_0   = vertex.get_band();
        //      int s_0   = vertex.get_e_spin();

        //      int b_1   = vertex.get_band();
        //      int s_1   = vertex.get_e_spin();

        int r_ind = vertex.get_r_site();
        int t_ind = t::dmn_size()/2+(vertex.get_tau()/delta_tau);

        if(t_ind<0 or t_ind>=t::dmn_size())
        t_ind = 0;

        double f_tau = 0.;
        if(vertex.get_HS_spin() != HS_ZERO)
        f_tau = 1.;

        K_r_t(0, 0, 0, 0, r_ind, t_ind) += sign*f_tau;
        K_r_t(0, 1, 0, 1, r_ind, t_ind) += sign*f_tau;
        }
      */
    }

    template<class parameters_type, class MOMS_type>
    int MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::find_first_non_interacting_spin(std::vector<vertex_singleton_type>& configuration_e_spin)
    {
      int configuration_size = configuration_e_spin.size();

      int vertex_index=0;
      while(vertex_index<configuration_size && configuration_e_spin[vertex_index].get_HS_spin() != HS_ZERO)
        vertex_index++;

      assert(vertex_index==configuration_size || configuration_e_spin[vertex_index].get_HS_spin() == HS_ZERO);

      return vertex_index;
    }

    template<class parameters_type, class MOMS_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::compute_M_r_w(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w)
    {
      function<std::complex<double>, dmn_2<w_dmn_t, p_dmn_t> > tmp("tmp M_r_w");
      cached_nfft_1D_M_r_w_obj.finalize(tmp);

      for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++)
          for(int s2_ind=0; s2_ind<s_dmn_t::dmn_size(); s2_ind++)
            for(int b2_ind=0; b2_ind<b_dmn_t::dmn_size(); b2_ind++)
              for(int s1_ind=0; s1_ind<s_dmn_t::dmn_size(); s1_ind++)
                for(int b1_ind=0; b1_ind<b_dmn_t::dmn_size(); b1_ind++)
                  M_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind)
                    = tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

      //double one_div_n_sites = 1./double(base_cluster_type::get_cluster_size());
      double one_div_n_sites = 1./double(r_DCA::dmn_size());
      M_r_w *= one_div_n_sites;
    }

    template<class parameters_type, class MOMS_type>
    void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::finalize(function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w,
                                                                                                   function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& M_r_w_squared)
    {
      {
        function<std::complex<double>, dmn_2<w_dmn_t, p_dmn_t> > tmp("tmp M_r_w");
        cached_nfft_1D_M_r_w_obj.finalize(tmp);

        for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
          for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++)
            for(int s2_ind=0; s2_ind<s_dmn_t::dmn_size(); s2_ind++)
              for(int b2_ind=0; b2_ind<b_dmn_t::dmn_size(); b2_ind++)
                for(int s1_ind=0; s1_ind<s_dmn_t::dmn_size(); s1_ind++)
                  for(int b1_ind=0; b1_ind<b_dmn_t::dmn_size(); b1_ind++)
                    M_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind)
                      = tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

        //double one_div_n_sites = 1./double(base_cluster_type::get_cluster_size());
        double one_div_n_sites = 1./double(r_DCA::dmn_size());
        M_r_w *= one_div_n_sites;
      }

      {
        function<std::complex<double>, dmn_2<w_dmn_t, p_dmn_t> > tmp("tmp M_r_w");
        cached_nfft_1D_M_r_w_squared_obj.finalize(tmp);

        for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
          for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++)
            for(int s2_ind=0; s2_ind<s_dmn_t::dmn_size(); s2_ind++)
              for(int b2_ind=0; b2_ind<b_dmn_t::dmn_size(); b2_ind++)
                for(int s1_ind=0; s1_ind<s_dmn_t::dmn_size(); s1_ind++)
                  for(int b1_ind=0; b1_ind<b_dmn_t::dmn_size(); b1_ind++)
                    M_r_w_squared(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind)
                      = tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

        //double one_div_n_sites = 1./double(base_cluster_type::get_cluster_size());
        double one_div_n_sites = 1./double(r_DCA::dmn_size());
        M_r_w_squared *= one_div_n_sites;
      }
    }

#ifdef MEASURE_ERROR_BARS
    //     /*!
    //      *  \brief Output and store standard deviation and error.
    //      *
    //      *  It computes and write to the given files the standard deviation of the measurements of the one particle accumulator.
    //      *  It outputs the L1-Norm, i.e. \f$\sum_{i=1}^N \left|x_i\right|/N\f$, the L2-Norm, i.e. \f$\sqrt{\sum_{i=1}^N \left|x_i\right|^2/N}\f$,
    //      *  and the Linf-Norm, i.e. \f$\max_{i=1}^N \left|x_i\right|\f$ of the standard deviation and of the error.
    //      */
    //     template<class parameters_type, class MOMS_type>
    //     void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::store_standard_deviation(int nr_measurements,
    //                                                                                                                    std::ofstream& points_file,
    //                                                                                                                    std::ofstream& norm_file)
    //     {
    //       std::pair<std::vector<scalar_type>, int> std = cached_nfft_1D_obj.get_standard_deviation(concurrency, nr_measurements, points_file, norm_file);
    //       scalar_type sqrt_n = sqrt(std.second);

    //       cout<<scientific;
    //       cout.precision(6);
    //       if(concurrency.id()==concurrency.first()){
    //         cout << "\n\n";
    //         cout << "\t\t             standard deviation ||  error \n";
    //         cout << "\t\t L1-norm   : " << std.first[0] << "       ||  " << std.first[0]/sqrt_n << "\n";
    //         cout << "\t\t L2-norm   : " << std.first[1] << "       ||  " << std.first[1]/sqrt_n << "\n";
    //         cout << "\t\t Linf-norm : " << std.first[2] << "       ||  " << std.first[2]/sqrt_n << "\n";
    //         cout << "\n\n";
    //       }
    //     }


    //     /*!
    //      *  \brief Update the sum of the squares of the measurements. It has to be called after each measurement.
    //      */
    //     template<class parameters_type, class MOMS_type>
    //     void MC_single_particle_accumulator<CT_AUX_SOLVER, NFFT, parameters_type, MOMS_type>::update_sum_squares()
    //     {
    //       cached_nfft_1D_obj.update_sum_squares();
    //     }
#endif

  }

}

#endif
