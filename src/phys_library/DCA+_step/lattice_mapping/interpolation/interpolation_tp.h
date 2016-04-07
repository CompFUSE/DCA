//-*-C++-*-

#ifndef DCA_INTERPOLATION_TP_H
#define DCA_INTERPOLATION_TP_H
#include"phys_library/domain_types.hpp"
#include "math_library/functional_transforms/basis_functions/basis_functions.hpp"
using namespace types;

namespace DCA
{
  /*!
   *  \author Peter Staar
   *  \brief  This class computes the interpolated cluster vertex.
   */
  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  class interpolation_tp : public interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>
  {

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef typename source_k_dmn::parameter_type::dual_type source_r_cluster_type;

    typedef dmn_0<centered_cluster_domain<source_r_cluster_type> >  r_centered_dmn;

    typedef math_algorithms::functional_transforms::basis_function<typename target_k_dmn::parameter_type, math_algorithms::KRONECKER_DELTA,
					    typename source_k_dmn::parameter_type, math_algorithms::HERMITE_CUBIC_SPLINE> basis_function_type;

  public:

    interpolation_tp(parameters_type& parameters_ref);
    ~interpolation_tp();

    template<typename scalartype>
    void initialize_T_K_to_k(LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T_K_to_k);

    template<typename scalartype>
    void execute(FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b,b,k_DCA        ,w_VERTEX>, dmn_4<b,b,k_DCA        ,w_VERTEX> > >& Gamma_cluster,
		 FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b,b,k_HOST_VERTEX,w_VERTEX>, dmn_4<b,b,k_HOST_VERTEX,w_VERTEX> > >& Gamma_lattice);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;
  };

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::interpolation_tp(parameters_type& parameters_ref):
    interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>(parameters_ref),

    parameters(parameters_ref),
    concurrency(parameters.get_concurrency())
  {
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::~interpolation_tp()
  {}

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  template<typename scalartype>
  void interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::initialize_T_K_to_k(LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T_K_to_k)
  {
    int Nr = target_k_dmn::dmn_size();
    int Nc = source_k_dmn::dmn_size();

    T_K_to_k.resize_no_copy(std::pair<int,int>(Nr,Nc));

    for(int j=0; j<Nc; j++)
      for(int i=0; i<Nr; i++)
	T_K_to_k(i,j) = basis_function_type::execute(i,j);
    
  }

  template<typename parameters_type, typename source_k_dmn, typename target_k_dmn>
  template<typename scalartype>
  void interpolation_tp<parameters_type, source_k_dmn, target_k_dmn>::execute(FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b,b,k_DCA        ,w_VERTEX>, dmn_4<b,b,k_DCA        ,w_VERTEX> > >& Gamma_cluster,
									      FUNC_LIB::function<std::complex<scalartype>, dmn_2<dmn_4<b,b,k_HOST_VERTEX,w_VERTEX>, dmn_4<b,b,k_HOST_VERTEX,w_VERTEX> > >& Gamma_lattice)
  {
    LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU> T_K_to_k("T_K_to_k");

    initialize_T_K_to_k(T_K_to_k);

    math_algorithms::functional_transforms::TRANSFORM<k_DCA, k_HOST_VERTEX>::execute_on_all(Gamma_cluster, Gamma_lattice, T_K_to_k);
  }
    
}

#endif
