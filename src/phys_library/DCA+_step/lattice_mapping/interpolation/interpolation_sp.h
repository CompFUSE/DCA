//-*-C++-*-

#ifndef DCA_INTERPOLATION_SP_H
#define DCA_INTERPOLATION_SP_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  /*!
   *  \author Peter Staar
   *  \brief  This class computes the interpolated cluster self-energy, using the alpha transformation.
   */
  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  class interpolation_sp : public interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>
  {

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef typename source_k_dmn::parameter_type::dual_type source_r_cluster_type;

    typedef dmn_0<centered_cluster_domain<source_r_cluster_type> >  r_centered_dmn;

    typedef dmn_3<nu, nu, r_centered_dmn>    nu_nu_r_centered;
    typedef dmn_4<nu, nu, r_centered_dmn, w> nu_nu_r_centered_w;

  public:

    interpolation_sp(parameters_type& parameters_ref);
    ~interpolation_sp();

    void execute_with_alpha_transformation(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,source_k_dmn,w> >& cluster_self_energy,
                                           FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,target_k_dmn,w> >& interp_self_energy);
  private:

    void execute(FUNC_LIB::function<std::complex<double>, dmn_3<nu,nu,source_k_dmn> >& cluster_self_energy,
                 FUNC_LIB::function<std::complex<double>, dmn_3<nu,nu,target_k_dmn> >& interp_self_energy);

    void execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,source_k_dmn,w> >& cluster_self_energy,
                 FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,target_k_dmn,w> >& interp_self_energy);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;
  };

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::interpolation_sp(parameters_type& parameters_ref):
    interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>(parameters_ref),

    parameters(parameters_ref),
    concurrency(parameters.get_concurrency())
  {
    if(concurrency.id() == concurrency.first())
      std::cout << "\n\n\t" << __FUNCTION__ << " is created " << print_time();
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::~interpolation_sp()
  {}

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute_with_alpha_transformation(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,source_k_dmn,w> >& cluster_self_energy,
                                                                                                        FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,target_k_dmn,w> >& interp_self_energy)
  {
    r_centered_dmn::parameter_type::initialize();

    FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w>  cluster_alpha_k     ("cluster_alpha_k");
    FUNC_LIB::function<std::complex<double>, nu_nu_k_HOST_w> interpolated_alpha_k("interpolated_alpha");

    transform_to_alpha::forward(1., cluster_self_energy, cluster_alpha_k);

    execute(cluster_alpha_k, interpolated_alpha_k);

    transform_to_alpha::backward(1., interp_self_energy, interpolated_alpha_k);
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute(FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,source_k_dmn,w> >& cluster_function,
                                                                              FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,target_k_dmn,w> >& interp_function)
  {
    r_centered_dmn::parameter_type::initialize();

    FUNC_LIB::function<std::complex<double>, nu_nu_r_centered_w> cluster_centered_function("cluster_centered_function");

    math_algorithms::functional_transforms::TRANSFORM<source_k_dmn, r_centered_dmn>::execute(cluster_function, cluster_centered_function);

    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
      for(int r_ind=0; r_ind<r_centered_dmn::dmn_size(); r_ind++)
        for(int j=0; j<nu::dmn_size(); j++)
          for(int i=0; i<nu::dmn_size(); i++)
            cluster_centered_function(i,j,r_ind,w_ind) *= r_centered_dmn::parameter_type::get_weights()[r_ind];

    math_algorithms::functional_transforms::TRANSFORM<r_centered_dmn, target_k_dmn>::execute(cluster_centered_function, interp_function);
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute(FUNC_LIB::function<std::complex<double>, dmn_3<nu,nu,source_k_dmn> >& cluster_function,
                                                                              FUNC_LIB::function<std::complex<double>, dmn_3<nu,nu,target_k_dmn> >& interp_function)
  {
    r_centered_dmn::parameter_type::initialize();

    FUNC_LIB::function<std::complex<double>, nu_nu_r_centered> cluster_centered_function("cluster_centered_function");

    math_algorithms::functional_transforms::TRANSFORM<source_k_dmn, r_centered_dmn>::execute(cluster_function, cluster_centered_function);

    for(int r_ind=0; r_ind<r_centered_dmn::dmn_size(); r_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          cluster_centered_function(i,j,r_ind) *= r_centered_dmn::parameter_type::get_weights()[r_ind];

    math_algorithms::functional_transforms::TRANSFORM<r_centered_dmn, target_k_dmn>::execute(cluster_centered_function, interp_function);
  }
  
}

#endif
