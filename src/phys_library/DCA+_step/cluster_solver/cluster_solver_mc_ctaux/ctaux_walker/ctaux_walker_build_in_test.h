//-*-C++-*-

#ifndef DCA_QMCI_CT_AUX_WALKER_BIT_H
#define DCA_QMCI_CT_AUX_WALKER_BIT_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  namespace QMCI
  {
    /*!
     *  \defgroup CT-AUX-WALKER
     *  \ingroup  CT-AUX
     */

    /*!
     *  \defgroup STRUCTURES
     *  \ingroup  CT-AUX
     */

    /*!
     *  \ingroup CT-AUX
     *
     *  \brief   This class organizes the MC-walker in the CT-AUX QMC
     *  \author  Peter Staar
     *  \version 1.0
     */
    template<class parameters_type, class MOMS_type>
    class MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>
    {

      typedef vertex_singleton                         vertex_singleton_type;
      typedef CT_AUX_HS_configuration<parameters_type> configuration_type;

      typedef typename parameters_type::random_number_generator rng_type;

      typedef typename MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type>::profiler_type    profiler_type;
      typedef typename MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type>::concurrency_type concurrency_type;

    public:

      MC_walker_BIT(parameters_type& parameters_ref,
                    MOMS_type&       MOMS_ref,
                    int              id);


      void  initialize();

      FUNC_LIB::function<double, dmn_0<numerical_error_domain> >& get_error_distribution();

      template<LIN_ALG::device_type device_t>
      void check_G0_matrices(configuration_type& configuration,
                             LIN_ALG::matrix<double, device_t>& G0_up,
                             LIN_ALG::matrix<double, device_t>& G0_dn);

      template<LIN_ALG::device_type device_t>
      void check_N_matrices(configuration_type& configuration,
                            LIN_ALG::matrix<double, device_t>& G0_up,
                            LIN_ALG::matrix<double, device_t>& G0_dn,
                            LIN_ALG::matrix<double, device_t>& N_up,
                            LIN_ALG::matrix<double, device_t>& N_dn);

      template<LIN_ALG::device_type device_t>
      void check_G_matrices(configuration_type& configuration,
                            LIN_ALG::matrix<double, device_t>& G0_up,
                            LIN_ALG::matrix<double, device_t>& G0_dn,
                            LIN_ALG::matrix<double, device_t>& N_up,
                            LIN_ALG::matrix<double, device_t>& N_dn,
                            LIN_ALG::matrix<double, device_t>& G_up,
                            LIN_ALG::matrix<double, device_t>& G_dn);

    private:

      parameters_type&   parameters;
      MOMS_type&         MOMS;
      concurrency_type&  concurrency;

      int                thread_id;

      CV<parameters_type> CV_obj;

      G0_INTERPOLATION<LIN_ALG::CPU, parameters_type> G0_CPU_tools_obj;
      N_TOOLS         <LIN_ALG::CPU, parameters_type> N_CPU_tools_obj;
      G_TOOLS         <LIN_ALG::CPU, parameters_type> G_CPU_tools_obj;

      LIN_ALG::matrix<double, LIN_ALG::CPU> G0_up_CPU;
      LIN_ALG::matrix<double, LIN_ALG::CPU> G0_dn_CPU;

      LIN_ALG::matrix<double, LIN_ALG::CPU> N_up_CPU;
      LIN_ALG::matrix<double, LIN_ALG::CPU> N_dn_CPU;

      LIN_ALG::matrix<double, LIN_ALG::CPU> G_up_CPU;
      LIN_ALG::matrix<double, LIN_ALG::CPU> G_dn_CPU;

      FUNC_LIB::function<double, dmn_0<numerical_error_domain> > error;
    };

    template<class parameters_type, class MOMS_type>
    MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::MC_walker_BIT(parameters_type& parameters_ref,
                                                                            MOMS_type&       MOMS_ref,
                                                                            int              id):
      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      CV_obj(parameters),

      G0_CPU_tools_obj(thread_id, parameters),
      N_CPU_tools_obj (thread_id, parameters, CV_obj),
      G_CPU_tools_obj (thread_id, parameters, CV_obj),

      G0_up_CPU("G0_up_CPU (MC_walker_BIT)"),
      G0_dn_CPU("G0_up_CPU (MC_walker_BIT)"),

      N_up_CPU("N_up_CPU (MC_walker_BIT)"),
      N_dn_CPU("N_up_CPU (MC_walker_BIT)"),

      G_up_CPU("G_up_CPU (MC_walker_BIT)"),
      G_dn_CPU("G_up_CPU (MC_walker_BIT)"),

      error("error")
    {}

    template<class parameters_type, class MOMS_type>
    FUNC_LIB::function<double, dmn_0<numerical_error_domain> >& MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::get_error_distribution()
    {
      return error;
    }

    template<class parameters_type, class MOMS_type>
    void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::initialize()
    {
      error = 0.;

      CV_obj          .initialize(MOMS);
      G0_CPU_tools_obj.initialize(MOMS);
    }

    template<class parameters_type, class MOMS_type>
    template<LIN_ALG::device_type device_t>
    void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_G0_matrices(configuration_type& configuration,
                                                                                     LIN_ALG::matrix<double, device_t>& G0_up,
                                                                                     LIN_ALG::matrix<double, device_t>& G0_dn)
    {
      G0_CPU_tools_obj.build_G0_matrix(configuration, G0_up_CPU, e_UP);
      G0_CPU_tools_obj.build_G0_matrix(configuration, G0_dn_CPU, e_DN);

      G0_up_CPU.difference(G0_up);
      G0_dn_CPU.difference(G0_dn);
    }

    template<class parameters_type, class MOMS_type>
    template<LIN_ALG::device_type device_t>
    void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_N_matrices(configuration_type& configuration,
                                                                                    LIN_ALG::matrix<double, device_t>& G0_up,
                                                                                    LIN_ALG::matrix<double, device_t>& G0_dn,
                                                                                    LIN_ALG::matrix<double, device_t>& N_up,
                                                                                    LIN_ALG::matrix<double, device_t>& N_dn)
    {
      G0_up_CPU.difference(G0_up);
      G0_dn_CPU.difference(G0_dn);

      N_CPU_tools_obj.build_N_matrix(configuration, N_up_CPU, G0_up_CPU, e_UP);
      N_CPU_tools_obj.build_N_matrix(configuration, N_dn_CPU, G0_dn_CPU, e_DN);

      double err_up = N_up_CPU.difference(N_up);
      double err_dn = N_dn_CPU.difference(N_dn);

      std::vector<double>& x = numerical_error_domain::get_elements();
      for(size_t l=0; l<x.size()-1; l++)
        if((err_up>x[l] and err_up<x[l+1]) or
           (err_dn>x[l] and err_dn<x[l+1]))
          error(l) += 1;
    }

    template<class parameters_type, class MOMS_type>
    template<LIN_ALG::device_type device_t>
    void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_G_matrices(configuration_type& configuration,
                                                                                    LIN_ALG::matrix<double, device_t>& G0_up,
                                                                                    LIN_ALG::matrix<double, device_t>& G0_dn,
                                                                                    LIN_ALG::matrix<double, device_t>& N_up,
                                                                                    LIN_ALG::matrix<double, device_t>& N_dn,
                                                                                    LIN_ALG::matrix<double, device_t>& G_up,
                                                                                    LIN_ALG::matrix<double, device_t>& G_dn)
    {
      G0_up_CPU.difference(G0_up);
      G0_dn_CPU.difference(G0_dn);

      N_up_CPU.difference(N_up);
      N_dn_CPU.difference(N_dn);

      G_CPU_tools_obj.build_G_matrix(configuration, N_up_CPU, G0_up_CPU, G_up_CPU, e_UP);
      G_CPU_tools_obj.build_G_matrix(configuration, N_dn_CPU, G0_dn_CPU, G_dn_CPU, e_DN);

      G_up_CPU.difference(G_up);
      G_dn_CPU.difference(G_dn);
    }

  }

}

#endif
