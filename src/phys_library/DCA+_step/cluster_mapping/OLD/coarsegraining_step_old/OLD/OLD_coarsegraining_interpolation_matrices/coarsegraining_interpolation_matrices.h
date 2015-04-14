//-*-C++-*-

#ifndef DCA_COARSEGRAIN_INTERPOLATION_MATRICES_H
#define DCA_COARSEGRAIN_INTERPOLATION_MATRICES_H

namespace DCA
{
  template<typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
  class interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME> > >
  {  
    typedef dmn_0<coarsegraining_domain<K_dmn, NAME> > q_dmn;
    typedef typename k_dmn::parameter_type::dual_type  r_dmn;

    typedef MATH_ALGORITHMS::basis_transform<typename k_dmn::parameter_type, r_dmn> trafo_k_to_r_type;
    typedef MATH_ALGORITHMS::basis_transform<r_dmn, typename q_dmn::parameter_type> trafo_r_to_q_type;

    typedef typename trafo_k_to_r_type::matrix_type trafo_matrix_type;

  public:

    typedef LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> matrix_type;

  public:
    
    static FUNC_LIB::function<matrix_type, K_dmn>& get();

    static matrix_type& get(int k_ind);

    static bool is_initialized();

    template<typename parameters_type>
    static void initialize(parameters_type& parameters);
  };

}

#endif
