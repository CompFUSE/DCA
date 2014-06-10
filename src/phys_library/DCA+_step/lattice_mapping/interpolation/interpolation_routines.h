//-*-C++-*-

#ifndef DCA_INTERPOLATION_ROUTINES_H
#define DCA_INTERPOLATION_ROUTINES_H

namespace DCA
{
  /*!
   *  \author Peter Staar
   *  \brief  This class computes the interpolated cluster self-energy, using the alpha transformation.
   */
  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  class interpolation_routines
  {
#include "type_definitions.h"

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef typename source_k_dmn::parameter_type::dual_type source_r_cluster_type;

    typedef dmn_0<centered_cluster_domain<source_r_cluster_type> >  r_centered_dmn;

    typedef dmn_3<nu, nu, r_centered_dmn>    nu_nu_r_centered;
    typedef dmn_4<nu, nu, r_centered_dmn, w> nu_nu_r_centered_w;

  public:

    interpolation_routines(parameters_type& parameters_ref);
    ~interpolation_routines();

  private:

    void initialize();

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    //MATH_LIBRARY::gaussian_fit<double, source_k_dmn, target_k_dmn> gaussian_fit_obj;
  };

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>::interpolation_routines(parameters_type& parameters_ref):
    parameters(parameters_ref),
    concurrency(parameters.get_concurrency())
  {
    initialize();

//     if(concurrency.id() == concurrency.first())
//       cout << "\n\n\t" << __FUNCTION__ << " is created " << print_time();
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>::~interpolation_routines()
  {}

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  void interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>::initialize()
  {
    //gaussian_fit_obj.initialize_K_to_k(true, 1.e-3);
  }

}

#endif
