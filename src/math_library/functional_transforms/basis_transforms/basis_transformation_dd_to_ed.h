//-*-C++-*-

#ifndef BASIS_TRANSFORMATION_DD_TO_ED_H
#define BASIS_TRANSFORMATION_DD_TO_ED_H

namespace MATH_ALGORITHMS
{
  /*! 
   *  \ingroup 
   *
   *  \brief  
   *  \author Peter Staar
   *
   *  \version 1.0
   *  \date    2013
   */
  template<typename input_type, typename output_type>
  class basis_transformation<input_type, DISCRETE, output_type, EXPANSION>
  {
  public:

    typedef input_type  rh_dmn_type;
    typedef output_type lh_dmn_type;

    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;

    typedef basis_transformation<output_type, EXPANSION, input_type, DISCRETE> inverse_basis_transformation_type;

    typedef typename inverse_basis_transformation_type::matrix_type matrix_type;

    //     typedef typename lh_spec_dmn_type::scalar_type f_scalar_type;
    //     typedef LIN_ALG::matrix<f_scalar_type, LIN_ALG::CPU> matrix_type;

  public:

    static bool& is_initialized()
    {
      static bool initialized = false;
      return initialized;
    }
    
    static std::string& get_name()
    {
      static std::string name = "basis-transformation";
      return name;
    }

    static matrix_type& get_transformation_matrix()
    {
      static matrix_type T;

      if(not is_initialized())
	initialize_transformation_matrix();

      return T;
    }

    static void initialize_transformation_matrix()
    {
      GENERIC_ASSERT<DISCRETE  == rh_spec_dmn_type::DOMAIN_REPRESENTATION>::execute();
      GENERIC_ASSERT<EXPANSION == lh_spec_dmn_type::DOMAIN_REPRESENTATION>::execute();

      is_initialized() = true;

      int M = lh_dmn_type::get_size();
      int N = rh_dmn_type::get_size();

      assert(M>0 and N>0);
 
      matrix_type& T = get_transformation_matrix();
      
      T.resize_no_copy(std::pair<int, int>(M, N));

      inverse_basis_transformation_type::is_initialized() = false;

      matrix_type& T_inv = inverse_basis_transformation_type::get_transformation_matrix();
      
      LIN_ALG::PSEUDO_INVERSE<LIN_ALG::CPU>::execute(T_inv, T);
    }

  };

}

#endif
