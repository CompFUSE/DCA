//-*-C++-*-

#ifndef COMPUTE_THIRD_ORDER_SIGMA_H
#define COMPUTE_THIRD_ORDER_SIGMA_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{

  namespace SERIES_EXPANSION
  {
    /*!
     * \authors Peter Staar, Urs R. Haehner
     *
     * \brief  This class implements the computation of the self-energy in second order.
     */
    template<class parameter_type, class k_dmn_t>
    class sigma_perturbation<3, parameter_type, k_dmn_t>
    {

    public:

      typedef typename compute_interaction::function_type U_type;

      typedef compute_bubble<ph, parameter_type, k_dmn_t, w> ph_bubble_t;
      typedef compute_bubble<ph, parameter_type, k_dmn_t, w> pp_bubble_t;

      typedef typename ph_bubble_t::function_type chi_type;
      typedef typename pp_bubble_t::function_type phi_type;

      typedef FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > function_type;

    public:

      sigma_perturbation(parameter_type&      parameters_ref,
                         compute_interaction& interaction_obj,
                         compute_bubble<ph, parameter_type, k_dmn_t, w>& chi_obj,
                         compute_bubble<pp, parameter_type, k_dmn_t, w>& phi_obj);

      ~sigma_perturbation();

      function_type& get_function() { return Sigma; }

      void execute_on_cluster(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G);

      template<IO::FORMAT DATA_FORMAT>
      void write(IO::writer<DATA_FORMAT>& writer);

    private:

      void execute_RPA(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G);
      void execute_VC(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G);

      int subtract_freq_fb(int, int);

    protected:

      parameter_type& parameters;

      U_type& U;

      chi_type& chi;
      phi_type& phi;

      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > Sigma;
      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > Sigma_RPA;
      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> > Sigma_VC;
    };

    template<class parameter_type, class k_dmn_t>
    sigma_perturbation<3, parameter_type, k_dmn_t>::sigma_perturbation(parameter_type&      parameters_ref,
                                                                       compute_interaction& interaction_obj,
                                                                       compute_bubble<ph, parameter_type, k_dmn_t, w>& chi_obj,
                                                                       compute_bubble<pp, parameter_type, k_dmn_t, w>& phi_obj):
      parameters(parameters_ref),

      U(interaction_obj.get_function()),

      chi(chi_obj.get_function()),
      phi(phi_obj.get_function()),

      Sigma("Sigma-3rd-order"),
      Sigma_RPA("Sigma-3rd-order-RPA"),
      Sigma_VC("Sigma-3rd-order-VC")
    {}

    template<class parameter_type, class k_dmn_t>
    sigma_perturbation<3, parameter_type, k_dmn_t>::~sigma_perturbation()
    {}

    template<class parameter_type, class k_dmn_t>
    template<IO::FORMAT DATA_FORMAT>
    void sigma_perturbation<3, parameter_type, k_dmn_t>::write(IO::writer<DATA_FORMAT>& /*writer*/)
    //WARNING empty function
    {

    }


    template<class parameter_type, class k_dmn_t>
    void sigma_perturbation<3, parameter_type, k_dmn_t>::execute_on_cluster(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G)
    {
      std::cout << __FUNCTION__ << std::endl;

      std::cout << "\t U : " << U(0,0,0,1) << std::endl;

      sigma_perturbation<3, parameter_type, k_dmn_t>::execute_RPA(G);
      sigma_perturbation<3, parameter_type, k_dmn_t>::execute_VC(G);

      Sigma = 0.;
      Sigma += Sigma_RPA;
      Sigma += Sigma_VC;
    }

    template<class parameter_type, class k_dmn_t>
    void sigma_perturbation<3, parameter_type, k_dmn_t>::execute_RPA(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G)
    {
      std::cout << __FUNCTION__ << std::endl;

      double U_value = U(0,0,0,1);

      Sigma_RPA = 0.;

      for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){
        for(int q_ind=0; q_ind<k_dmn_t::dmn_size(); ++q_ind){

          int nu_c = (nu_ind-w_VERTEX_BOSONIC::dmn_size()/2);

          for(int w_ind=std::fabs(nu_c); w_ind<w::dmn_size()-std::fabs(nu_c); ++w_ind){
            for (int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){
              int k_minus_q  = k_dmn_t::parameter_type::subtract(q_ind, k_ind);
              int w_minus_nu = w_ind-nu_c;
              Sigma_RPA(0,0, 0,0, k_ind, w_ind) += G(0,0, 0,0, k_minus_q, w_minus_nu) * chi(0,0, 0,0, q_ind, nu_ind) * chi(0,0, 0,0, q_ind, nu_ind);
            }
          }
        }
      }
      for(int w_ind=0; w_ind<w::dmn_size(); ++w_ind)
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind)
          Sigma_RPA(0,1, 0,1, k_ind, w_ind) = Sigma_RPA(0,0, 0,0, k_ind, w_ind);

      double factor = 1./(parameters.get_beta()*k_dmn_t::dmn_size())*U_value*U_value*U_value;
      Sigma_RPA *= factor;
    }

    template<class parameter_type, class k_dmn_t>
    void sigma_perturbation<3, parameter_type, k_dmn_t>::execute_VC(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_dmn_t, w> >& G)
    {
      std::cout << __FUNCTION__ << std::endl;

      double U_value = U(0,0,0,1);

      Sigma_VC = 0.;

      for(int nu_ind=0; nu_ind<w_VERTEX_BOSONIC::dmn_size(); ++nu_ind){
        for(int q_ind=0; q_ind<k_dmn_t::dmn_size(); ++q_ind){

          for(int w_ind=0; w_ind<w::dmn_size(); ++w_ind){
            int nu_minus_w = subtract_freq_fb(w_ind, nu_ind);
            if (nu_minus_w<0 || nu_minus_w >= w::dmn_size()) continue;

            for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind){
              int q_minus_k = k_dmn_t::parameter_type::subtract(k_ind, q_ind);

              Sigma_VC(0,0, 0,0, k_ind, w_ind) += G(0,0, 0,0, q_minus_k, nu_minus_w) * phi(0,0, 0,0, q_ind, nu_ind) * phi(0,0, 0,0, q_ind, nu_ind);
            }
          }
        }
      }

      for(int w_ind=0; w_ind<w::dmn_size(); ++w_ind)
        for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind)
          Sigma_VC(0,1, 0,1, k_ind, w_ind) = Sigma_VC(0,0, 0,0, k_ind, w_ind);

      double factor = 1./(parameters.get_beta()*k_dmn_t::dmn_size())*U_value*U_value*U_value;
      Sigma_VC *= factor;
    }

    template<class parameter_type, class k_dmn_t>
    int sigma_perturbation<3, parameter_type, k_dmn_t>::subtract_freq_fb(int w1, int w2)
    {
      int w_f = 2*(w1 - w::dmn_size()/2) + 1;             // transform fermionic
      int w_b = 2*(w2 - w_VERTEX_BOSONIC::dmn_size()/2);  // transform bosonic
      int res = ((w_b-w_f) - 1 + w::dmn_size()) / 2;      // result is fermionic
      return res;
    }

  }

}

#endif
