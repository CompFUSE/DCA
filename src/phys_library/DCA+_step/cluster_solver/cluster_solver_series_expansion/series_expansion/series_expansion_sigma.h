//-*-C++-*-

#ifndef COMPUTE_SERIES_EXPANSION_SIGMA_H
#define COMPUTE_SERIES_EXPANSION_SIGMA_H

namespace DCA
{

  namespace SERIES_EXPANSION
  {
    /*!
     * \class
     *
     * \authors Peter Staar, Urs R. Haehner
     *
     * \brief  This class implements the computation of the self-energy using a series expansion.
     */
    template<class parameters_type, class MOMS_type>
    class series_expansion
    {
#include "type_definitions.h"

      typedef typename parameters_type::concurrency_type concurrency_type;

//       typedef typename MOMS_type::r_dmn_t r_dmn_t;
//       typedef typename MOMS_type::k_dmn_t k_dmn_t;

      typedef r_HOST r_dmn_t;
      typedef k_HOST k_dmn_t;

      typedef FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > sigma_function_t;

    public:

      series_expansion(parameters_type& parameter_ref,
                       MOMS_type&       MOMS_ref);
      ~series_expansion();

      template<class stream_type>
      void to_JSON(stream_type& ss);

      void execute(bool do_not_adjust_mu=true);

      sigma_function_t& get_Sigma() { return Sigma; }

      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

    private:

      parameters_type&  parameters;
      concurrency_type& concurrency;
      MOMS_type&        MOMS;

      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > Sigma;

      compute_interaction interaction_obj;

      compute_bubble<ph, parameters_type, k_dmn_t, w> ph_bubble;
      compute_bubble<pp, parameters_type, k_dmn_t, w> pp_bubble;

      sigma_perturbation<1, parameters_type, k_dmn_t> sigma_perturbation_1_obj;
      sigma_perturbation<2, parameters_type, k_dmn_t> sigma_perturbation_2_obj;
      sigma_perturbation<3, parameters_type, k_dmn_t> sigma_perturbation_3_obj;
      sigma_perturbation<4, parameters_type, k_dmn_t> sigma_perturbation_4_obj;
    };

    template<class parameters_type, class MOMS_type>
    series_expansion<parameters_type, MOMS_type>::series_expansion(parameters_type& parameters_ref,
                                                                   MOMS_type&       MOMS_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      MOMS(MOMS_ref),

      Sigma("perturbation-Sigma"),

      interaction_obj(),

      ph_bubble(parameters),
      pp_bubble(parameters),

      sigma_perturbation_1_obj(parameters, interaction_obj),
      sigma_perturbation_2_obj(parameters, interaction_obj, ph_bubble, pp_bubble),
      sigma_perturbation_3_obj(parameters, interaction_obj, ph_bubble, pp_bubble),
      sigma_perturbation_4_obj(parameters, interaction_obj, ph_bubble, pp_bubble)
    {
      interaction_obj.execute(MOMS.H_interactions);
    }

    template<class parameters_type, class MOMS_type>
    series_expansion<parameters_type, MOMS_type>::~series_expansion()
    {}

    template<class parameters_type, class MOMS_type>
    void series_expansion<parameters_type, MOMS_type>::execute(bool do_not_adjust_mu)
    {
//       if(do_not_adjust_mu)
//         MOMS.adjust_chemical_potential();

      //FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G_k_w = MOMS.G_k_w;
      //FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G_k_w = MOMS.G0_k_w_cluster_excluded;

      compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w> compute_lattice_Greens_function_obj(parameters, MOMS);

      compute_lattice_Greens_function_obj.execute();

      //FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_HOST, w> >& G_k_w = compute_lattice_Greens_function_obj.get_G0_k_w();
      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_HOST, w> >& G_k_w = compute_lattice_Greens_function_obj.get_G_k_w();

      //ph_bubble.execute_on_cluster(G_k_w);
      //pp_bubble.execute_on_cluster(G_k_w);

      ph_bubble.threaded_execute_on_cluster(G_k_w);
      pp_bubble.threaded_execute_on_cluster(G_k_w);

      //sigma_perturbation_1_obj.execute_on_cluster(MOMS.orbital_occupancy);
      //sigma_perturbation_2_obj.execute_on_cluster(G_k_w);
      //sigma_perturbation_3_obj.execute_on_cluster(G_k_w);
      //sigma_perturbation_4_obj.execute_on_cluster(G_k_w, sigma_perturbation_2_obj.get_function());

      sigma_perturbation_1_obj.         execute_on_cluster(MOMS.orbital_occupancy);
      sigma_perturbation_2_obj.threaded_execute_on_cluster(G_k_w);

      Sigma = 0.;

      Sigma += sigma_perturbation_1_obj.get_function();
      Sigma += sigma_perturbation_2_obj.get_function();
      //Sigma += sigma_perturbation_3_obj.get_function();
      //Sigma += sigma_perturbation_4_obj.get_function();

      if(true)
        {
          std::complex<double> I(0,1);
          for(int b_ind=0; b_ind<2*b::dmn_size(); ++b_ind){

            for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){

              int wc_ind = w::dmn_size()/8;

              double wc = w::get_elements()[wc_ind];

              std::complex<double> Sigma_wc = Sigma(b_ind, b_ind, k_ind, wc_ind);

              double alpha = real(Sigma_wc);
              double beta  = imag(Sigma_wc*wc);

              for(int w_ind=0; w_ind<wc_ind; ++w_ind){
                Sigma(b_ind, b_ind, k_ind,                 w_ind) = alpha + beta*I/w::get_elements()[w_ind];
                Sigma(b_ind, b_ind, k_ind, w::dmn_size()-1-w_ind) = alpha - beta*I/w::get_elements()[w_ind];
              }
            }
          }
        }
    }

    template<class parameters_type, class MOMS_type>
    template<IO::FORMAT DATA_FORMAT>
    void series_expansion<parameters_type, MOMS_type>::write(IO::writer<DATA_FORMAT>& writer)
    {
      writer.execute(Sigma);

      ph_bubble.write(writer);
      pp_bubble.write(writer);

      sigma_perturbation_1_obj.write(writer);
      sigma_perturbation_2_obj.write(writer);
      sigma_perturbation_3_obj.write(writer);
      sigma_perturbation_4_obj.write(writer);
    }
    
  }

}

#endif
