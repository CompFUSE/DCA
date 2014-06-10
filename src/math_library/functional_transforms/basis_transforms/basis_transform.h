//-*-C++-*-

#ifndef BASIS_TRANSFORM_H
#define BASIS_TRANSFORM_H

namespace MATH_ALGORITHMS
{
  template<typename input_type, typename output_type>
  class basis_transform
  {
  public:

    typedef typename input_type ::dmn_specifications_type input_spec_dmn_type;
    typedef typename output_type::dmn_specifications_type output_spec_dmn_type;

    typedef basis_transformation<input_type , input_spec_dmn_type ::DOMAIN_REPRESENTATION, 
				 output_type, output_spec_dmn_type::DOMAIN_REPRESENTATION> basis_transformation_type;

    typedef typename basis_transformation_type::matrix_type matrix_type;

  public:

    static bool& is_initialized()
    {
      return basis_transformation_type::is_initialized();
    }

    static std::string& get_name()
    {
      return basis_transformation_type::get_name();
    }

    static matrix_type& get_transformation_matrix()
    {
      return basis_transformation_type::get_transformation_matrix();
    }

  };


  template<typename input_type, typename output_type>
  class basis_transform<dmn_0<input_type>, dmn_0<output_type> >
  {
  public:

    typedef typename input_type ::dmn_specifications_type input_spec_dmn_type;
    typedef typename output_type::dmn_specifications_type output_spec_dmn_type;

    typedef basis_transformation<input_type , input_spec_dmn_type ::DOMAIN_REPRESENTATION, 
				 output_type, output_spec_dmn_type::DOMAIN_REPRESENTATION> basis_transformation_type;

    typedef typename basis_transformation_type::matrix_type matrix_type;

  public:

    static bool& is_initialized()
    {
      return basis_transformation_type::is_initialized();
    }

    static std::string& get_name()
    {
      return basis_transformation_type::get_name();
    }

    static matrix_type& get_transformation_matrix()
    {
      return basis_transformation_type::get_transformation_matrix();
    }
  };

}

#endif
