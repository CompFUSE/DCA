//-*-C++-*-

#ifndef DCA_QMCI_G0_INTERPOLATION_TEMPLATE_H
#define DCA_QMCI_G0_INTERPOLATION_TEMPLATE_H
#include "phys_library/domain_types.hpp"
#include "math_library/interpolation_library/akima_interpolation.h"
using namespace types;

namespace DCA
{
  namespace QMCI
  {
//     template<LIN_ALG::device_type device_t, typename parameters_type>
//     class G0_INTERPOLATION
//     {};

    /*!
     *  \class   G0_INTERPOLATION_TEMPLATE
     *  \ingroup CT-AUX-WALKER
     *
     *  \author Peter Staar
     *  \brief  This class organizes the interpolation of \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.
     */
    template<typename parameters_type>
    class G0_INTERPOLATION_TEMPLATE
    {

      typedef vertex_singleton    vertex_singleton_type;

      typedef r_DCA r_dmn_t;
      typedef k_DCA k_dmn_t;

      typedef typename r_dmn_t::parameter_type r_cluster_type;
      typedef typename k_dmn_t::parameter_type k_cluster_type;

      typedef typename parameters_type::concurrency_type concurrency_type;
      typedef typename parameters_type::profiler_type    profiler_t;

      typedef dmn_0<time_domain_left_oriented>                shifted_t;
      typedef dmn_4<nu, nu, r_dmn_t, shifted_t> nu_nu_r_dmn_t_shifted_t;

      typedef dmn_0<dmn<4, int> >                            akima_dmn_t;
      typedef dmn_5<akima_dmn_t, nu, nu, r_dmn_t, shifted_t> akima_nu_nu_r_dmn_t_shifted_t;

    public:

      G0_INTERPOLATION_TEMPLATE(int              id,
                                parameters_type& parameters);

      ~G0_INTERPOLATION_TEMPLATE();

      template<class MOMS_type>
      void initialize(MOMS_type&  MOMS);

    protected:

      template<class MOMS_type>
      void initialize_linear_coefficients(MOMS_type&  MOMS);

      template<class MOMS_type>
      void initialize_akima_coefficients(MOMS_type&  MOMS);

    protected:

      int thread_id;

      parameters_type&  parameters;
      concurrency_type& concurrency;

      nu_nu_r_dmn_t_shifted_t nu_nu_r_dmn_t_t_shifted_dmn;

      LIN_ALG::matrix<double, LIN_ALG::CPU> r1_minus_r0;

      FUNC_LIB::function<double, nu_nu_r_dmn_t_shifted_t>      G0_r_t_shifted;
      FUNC_LIB::function<double, nu_nu_r_dmn_t_shifted_t> grad_G0_r_t_shifted;

      FUNC_LIB::function<double, akima_nu_nu_r_dmn_t_shifted_t> akima_coefficients;

      int    N_t, linind, t_ind;
      double beta, N_div_beta, new_tau, scaled_tau, delta_tau, f_0, grad;
    };

    template<typename parameters_type>
    G0_INTERPOLATION_TEMPLATE<parameters_type>::G0_INTERPOLATION_TEMPLATE(int              id,
                                                                          parameters_type& parameters_ref):
      thread_id(id),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      r1_minus_r0(r_dmn_t::dmn_size())
    {
      beta       = parameters.get_beta();

      //N_t        = parameters.get_number_of_positive_times()+1;
      N_t        = parameters.get_sp_time_intervals()+1;
      N_div_beta = parameters.get_sp_time_intervals()/beta;

      for(int r1_ind=0; r1_ind<r_dmn_t::dmn_size(); r1_ind++)
        for(int r0_ind=0; r0_ind<r_dmn_t::dmn_size(); r0_ind++)
          r1_minus_r0(r0_ind, r1_ind) = r_cluster_type::subtract(r0_ind, r1_ind);
    }

    template<typename parameters_type>
    G0_INTERPOLATION_TEMPLATE<parameters_type>::~G0_INTERPOLATION_TEMPLATE()
    {}

    /*!
     *  \brief  Set the functions 'G0_r_t_shifted' and 'grad_G0_r_t_shifted'
     */
    template<typename parameters_type>
    template<class MOMS_type>
    void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize(MOMS_type&  MOMS)
    {
      initialize_linear_coefficients(MOMS);

      initialize_akima_coefficients(MOMS);
    }

    template<typename parameters_type>
    template<class MOMS_type>
    void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize_linear_coefficients(MOMS_type&  MOMS)
    {
      for(int t_ind=0; t_ind<t::dmn_size()/2-1; t_ind++){

        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int nu1_ind=0; nu1_ind<b::dmn_size()*s::dmn_size(); nu1_ind++){
            for(int nu0_ind=0; nu0_ind<b::dmn_size()*s::dmn_size(); nu0_ind++){

              G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind) = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);
              grad_G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind) = (MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind+1)
                                                                     -MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind));
            }
          }
        }
      }

      for(int t_ind=t::dmn_size()/2; t_ind<t::dmn_size()-1; t_ind++){

        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int nu1_ind=0; nu1_ind<b::dmn_size()*s::dmn_size(); nu1_ind++){
            for(int nu0_ind=0; nu0_ind<b::dmn_size()*s::dmn_size(); nu0_ind++){
              G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind-1) = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);
              grad_G0_r_t_shifted(nu0_ind, nu1_ind, r_ind, t_ind-1) = (MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind+1)
                                                                       -MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind));
            }
          }
        }
      }
    }

    template<typename parameters_type>
    template<class MOMS_type>
    void G0_INTERPOLATION_TEMPLATE<parameters_type>::initialize_akima_coefficients(MOMS_type&  MOMS)
    {
      int size = t::dmn_size()/2;

      math_algorithms::interpolation::akima_interpolation<double> ai_obj(size);

      double* x = new double[size];
      double* y = new double[size];

      for(int t_ind=0; t_ind<t::dmn_size()/2; t_ind++)
        x[t_ind] = t_ind;

      {
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int nu1_ind=0; nu1_ind<b::dmn_size()*s::dmn_size(); nu1_ind++){
            for(int nu0_ind=0; nu0_ind<b::dmn_size()*s::dmn_size(); nu0_ind++){

              for(int t_ind=0; t_ind<t::dmn_size()/2; t_ind++)
                y[t_ind] = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

              ai_obj.initialize(x, y);

              for(int t_ind=0; t_ind<t::dmn_size()/2-1; t_ind++)
                for(int l=0; l<4; l++)
                  akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind) = ai_obj.get_alpha(l, t_ind);
            }
          }
        }
      }

      {
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int nu1_ind=0; nu1_ind<b::dmn_size()*s::dmn_size(); nu1_ind++){
            for(int nu0_ind=0; nu0_ind<b::dmn_size()*s::dmn_size(); nu0_ind++){

              for(int t_ind=t::dmn_size()/2; t_ind<t::dmn_size(); t_ind++)
                y[t_ind-t::dmn_size()/2] = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

              ai_obj.initialize(x, y);

              for(int t_ind=t::dmn_size()/2; t_ind<t::dmn_size()-1; t_ind++)
                for(int l=0; l<4; l++)
                  akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind-1) = ai_obj.get_alpha(l, t_ind-t::dmn_size()/2);
            }
          }
        }
      }

      delete [] x;
      delete [] y;

      /*
        {
        cout << "\n\n\n";
        for(int t_ind=0; t_ind<shifted_t::dmn_size(); t_ind++)
        cout << t_ind << "\t" << G0_r_t_shifted(0, 0, 0, t_ind) << endl;
        cout << "\n\n\n";

        for(int t_ind=0; t_ind<shifted_t::dmn_size(); t_ind++)
        {
        int linind    = 4*nu_nu_r_dmn_t_t_shifted_dmn(0,0,0,t_ind);
        double* a_ptr = &akima_coefficents(linind);

        for(double x=0; x<1.05; x+=0.1)
        cout << t_ind+x << "\t" << (a_ptr[0] + x*(a_ptr[1] + x*(a_ptr[2] + x*a_ptr[3]))) << endl;
        }


        sleep(10);
        throw std::logic_error(__FUNCTION__);
        }
      */
    }

  }

}

#endif
