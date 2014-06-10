//-*-C++-*-

#ifndef SS_CT_HYB_ACCUMULATOR_H
#define SS_CT_HYB_ACCUMULATOR_H

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \defgroup SS CT-HYB-ACCUMULATOR
     *  \ingroup  SS CT-HYB
     */

    /*!
     *  \ingroup SS CT-HYB
     *
     *  \brief   This class organizes the measurements in the SS CT-HYB QMC
     *  \author  Bart Ydens
     *  \version 1.0
     */
    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    class MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type> : public MC_accumulator_data,
                                                                            public ss_hybridization_solver_routines<parameters_type, MOMS_type>
    {
#include "type_definitions.h"

    public:

      typedef parameters_type   my_parameters_type;
      typedef MOMS_type         my_MOMS_type;

      typedef MC_walker<SS_CT_HYB, device_t, parameters_type, MOMS_type> walker_type;

      typedef ss_hybridization_solver_routines<parameters_type, MOMS_type> ss_hybridization_solver_routines_type;

      typedef typename walker_type::ss_hybridization_walker_routines_type ss_hybridization_walker_routines_type;

      typedef r_DCA r_dmn_t;
      typedef k_DCA k_dmn_t;

      typedef w                      w_dmn_t;
      typedef dmn_3<nu, nu, r_dmn_t> p_dmn_t;

      typedef b b_dmn_t;
      typedef s s_dmn_t;

      typedef typename parameters_type::profiler_type    profiler_type;
      typedef typename parameters_type::concurrency_type concurrency_type;

      typedef double scalar_type;

      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::vertex_vertex_matrix_type  vertex_vertex_matrix_type;
      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::orbital_configuration_type orbital_configuration_type;

      typedef typename MC_type_definitions<SS_CT_HYB, parameters_type, MOMS_type>::configuration_type configuration_type;

      typedef function<vertex_vertex_matrix_type, nu> M_matrix_type;

    public:

      MC_accumulator(parameters_type& parameters_ref,
                     MOMS_type&       MOMS_ref,
		     int              id=0);

      ~MC_accumulator();

      void initialize(int dca_iteration);

      void finalize();//function<double, nu> mu_DC);

      void measure(walker_type& walker);

      void update_from(walker_type& walker);
      void measure();

      void compute_G_r_w(function<double, nu> mu_DC);

      configuration_type& get_configuration() { return configuration; }

      function<double, dmn_0<Feynman_expansion_order_domain> >& get_visited_expansion_order_k() {return visited_expansion_order_k; }

      function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& get_G_r_w()  { return G_r_w; }
      function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> >& get_GS_r_w() { return GS_r_w; }

      void accumulate_length (walker_type&  walker);
      void accumulate_overlap(walker_type&  walker);

      function<double, nu>& get_length() { return length; }

      function<double, nu_nu>& get_overlap() { return overlap; }

      /*!
       *  \brief Print the functions G_r_w and G_k_w.
       */
      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

    protected:

      using MC_accumulator_data::DCA_iteration;
      using MC_accumulator_data::number_of_measurements;

      using MC_accumulator_data::current_sign;
      using MC_accumulator_data::accumulated_sign;

      parameters_type&  parameters;
      MOMS_type&        MOMS;
      concurrency_type& concurrency;

      int               thread_id;

      configuration_type                      configuration;
      function<vertex_vertex_matrix_type, nu> M_matrices;

      function<double, dmn_0<Feynman_expansion_order_domain> > visited_expansion_order_k;

      function<double ,       nu     > length;
      function<double , dmn_2<nu,nu> > overlap;

      function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> > G_r_w;
      function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w> > GS_r_w;

      MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, MOMS_type> single_particle_accumulator_obj;
    };

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::MC_accumulator(parameters_type& parameters_ref,
                                                                                    MOMS_type&       MOMS_ref,
										    int              id):
      ss_hybridization_solver_routines<parameters_type, MOMS_type>(parameters_ref, MOMS_ref),

      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      configuration(),
      M_matrices("accumulator-M-matrices"),

      visited_expansion_order_k("visited-expansion-order-k"),

      length("length"),
      overlap("overlap"),

      G_r_w("G-r-w-measured"),
      GS_r_w("GS-r-w-measured"),

      single_particle_accumulator_obj(parameters)
    {}

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::~MC_accumulator()
    {}

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::initialize(int dca_iteration)
    {
      /*
        {
        for(int b_i=0; b_i<b::dmn_size(); b_i++){
        is_interacting_band_vector[b_i] = false;

        for(int b_j=0; b_j<parameters.get_interacting_bands().size(); b_j++)
        if(parameters.get_interacting_bands()[b_j]==b_i)
        is_interacting_band_vector[b_i] = true;
        }

        for(int b_i=0; b_i<b::dmn_size(); b_i++){
        cout << "\t band " << b_i << " is " << (is_interacting_band_vector[b_i]? "interacting\n" : "non-interacting\n");
        }
        }
      */

      MC_accumulator_data::initialize(dca_iteration);

      visited_expansion_order_k = 0;

      single_particle_accumulator_obj.initialize(G_r_w, GS_r_w);

      length  = 0;
      overlap = 0;
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::finalize()//function<double, nu> mu_DC)
    {
      single_particle_accumulator_obj.finalize(G_r_w, GS_r_w);

      //single_particle_accumulator_obj.compute(G_r_w, GS_r_w);

      /*
      {
        concurrency.sum_and_average(G_r_w  , -parameters.get_beta());
        concurrency.sum_and_average(GS_r_w , -parameters.get_beta());

        concurrency.sum_and_average(length , parameters.get_beta()*parameters.get_number_of_measurements());
        concurrency.sum_and_average(overlap , parameters.get_beta()*parameters.get_number_of_measurements());
      }
      */

      /*
      {
        for(int b_ind=0; b_ind<b::dmn_size(); b_ind++)
          if(ss_hybridization_solver_routines_type::is_interacting_band(b_ind))
            for(int s_ind=0; s_ind<s::dmn_size(); s_ind++)
              for(int w=0; w<w::dmn_size(); w++)
                GS_r_w(b_ind, s_ind,b_ind, s_ind,0,w) = GS_r_w(b_ind, s_ind,b_ind, s_ind,0,w);// - mu_DC(b_ind, s_ind)*G_r_w(b_ind, s_ind,b_ind, s_ind,0,w);
      }
      */
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    template<IO::FORMAT DATA_FORMAT>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& writer)
    {
      writer.execute(G_r_w);
      writer.execute(GS_r_w);
    }

    /*************************************************************
     **                                                         **
     **                    G2 - MEASUREMENTS                    **
     **                                                         **
     *************************************************************/

    /*
    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::measure(walker_type& walker)
    {
      number_of_measurements += 1;
      accumulated_sign       += current_sign;

      int k = walker.get_configuration().size();
      if(k < visited_expansion_order_k.size())
        visited_expansion_order_k(k) += 1;

      single_particle_accumulator_obj.accumulate(walker, MOMS.H_interactions);

      accumulate_length(walker);
      // accumulate_overlap(walker);
    }
    */

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::update_from(walker_type& walker)
    {
      current_sign = walker.get_sign();
      
      configuration.copy_from(walker.get_configuration());
      
      for(int l=0; l<nu::dmn_size(); l++)
        M_matrices(l).copy_from(walker.get_M_matrices()(l));
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::measure()
    {
      number_of_measurements += 1;
      accumulated_sign       += current_sign;

      int k = configuration.size();
      if(k < visited_expansion_order_k.size())
        visited_expansion_order_k(k) += 1;

      single_particle_accumulator_obj.accumulate(current_sign, configuration, M_matrices, MOMS.H_interactions);

      //accumulate_length();
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::accumulate_length(walker_type&  walker)
    {
      ss_hybridization_walker_routines_type& hybridization_routines = walker.get_ss_hybridization_walker_routines();

      Hybridization_vertex full_segment(0,parameters.get_beta());

      for(int ind=0; ind<b::dmn_size()*s::dmn_size(); ind++){
        length(ind) += hybridization_routines.compute_overlap(full_segment,
                                                              walker.get_configuration().get_vertices(ind),
                                                              walker.get_configuration().get_full_line(ind),
                                                              parameters.get_beta());
      }
    }

    template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
    void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::accumulate_overlap(walker_type&  walker)
    {
      ss_hybridization_walker_routines_type& hybridization_routines = walker.get_ss_hybridization_walker_routines();

      Hybridization_vertex full_segment(0,parameters.get_beta());

      for(int ind_1=0; ind_1<b::dmn_size()*s::dmn_size(); ind_1++){
        for(int ind_2=0; ind_2<b::dmn_size()*s::dmn_size(); ind_2++){
          if(walker.get_configuration().get_full_line(ind_1)){
            overlap(ind_1,ind_2) += hybridization_routines.compute_overlap(full_segment,
                                                                           walker.get_configuration().get_vertices(ind_2),
                                                                           walker.get_configuration().get_full_line(ind_2),
                                                                           parameters.get_beta());
          }
          else{
            for (typename orbital_configuration_type::iterator it=walker.get_configuration().get_vertices(ind_1).begin();
                 it!=walker.get_configuration().get_vertices(ind_1).end(); it++) {
              overlap(ind_1,ind_2) += hybridization_routines.compute_overlap(*it,
                                                                             walker.get_configuration().get_vertices(ind_2),
                                                                             walker.get_configuration().get_full_line(ind_2),
                                                                             parameters.get_beta());
            }
          }
        }
      }
    }


    /*
      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::accumulate_length()
      {
      //ss_hybridization_walker_routines_type& hybridization_routines = walker.get_ss_hybridization_walker_routines();

      Hybridization_vertex full_segment(0,parameters.get_beta());

      for(int ind=0; ind<b::dmn_size()*s::dmn_size(); ind++){
      length(ind) += hybridization_routines.compute_overlap(full_segment,
      configuration.get_vertices(ind),
      configuration.get_full_line(ind),
      parameters.get_beta());
      }
      }

      template<LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
      void MC_accumulator<SS_CT_HYB, device_t, parameters_type, MOMS_type>::accumulate_overlap()
      {
      //ss_hybridization_walker_routines_type& hybridization_routines = walker.get_ss_hybridization_walker_routines();

      Hybridization_vertex full_segment(0,parameters.get_beta());

      for(int ind_1=0; ind_1<b::dmn_size()*s::dmn_size(); ind_1++){
      for(int ind_2=0; ind_2<b::dmn_size()*s::dmn_size(); ind_2++){
      if(walker.get_configuration().get_full_line(ind_1))
      {
      overlap(ind_1,ind_2) += hybridization_routines.compute_overlap(full_segment,
      configuration.get_vertices(ind_2),
      configuration.get_full_line(ind_2),
      parameters.get_beta());
      }
      else
      {
      for (typename orbital_configuration_type::iterator it=configuration.get_vertices(ind_1).begin();
      it!=configuration.get_vertices(ind_1).end(); it++) {
      overlap(ind_1,ind_2) += hybridization_routines.compute_overlap(*it,
      configuration.get_vertices(ind_2),
      configuration.get_full_line(ind_2),
      parameters.get_beta());
      }
      }
      }
      }
      }
    */
  }
}

#endif

