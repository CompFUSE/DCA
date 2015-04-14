//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_ED_TO_CD_H
#define BASIS_TRANSFORMATIONS_ED_TO_CD_H

namespace MATH_ALGORITHMS
{
  template<typename type_input, typename type_output, int DMN_INDEX>
  class TRANSFORM_DOMAIN<type_input, EXPANSION, type_output, CONTINUOUS, DMN_INDEX> : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>
  {
  private:

    const static bool VERBOSE = false;

    typedef basis_transformation<type_input, EXPANSION, type_output, CONTINUOUS> basis_transformation_type;
    typedef typename basis_transformation_type::matrix_type                      matrix_type;

  public:
    
    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
			FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {      
      default_execute(f_input, f_output);
    }

  private:
    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void default_execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
				FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {

      if(VERBOSE)
	cout << "\n\t default-transform (expansion -> continuous) " << DMN_INDEX << "  " << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

      matrix_type& T = basis_transformation_type::get_transformation_matrix();

      TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
    }
 
  };

}

#endif
