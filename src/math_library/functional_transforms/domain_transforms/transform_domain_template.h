//-*-C++-*-

#ifndef DOMAIN_TRANSFORMATION_TEMPLATE_H
#define DOMAIN_TRANSFORMATION_TEMPLATE_H

namespace MATH_ALGORITHMS
{
  template<typename type_input,  DOMAIN_REPRESENTATIONS DMN_REP_LHS,
	   typename type_output, DOMAIN_REPRESENTATIONS DMN_REP_RHS,
	   int DMN_INDEX>
  struct TRANSFORM_DOMAIN
  {
    const static bool VERBOSE = true;
    
    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
			FUNC_LIB::function<scalartype_output, domain_output>& f_output);
  };

}

#endif
