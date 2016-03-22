//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_DMN_2_TO_Y_LM_DMN_H
#define BASIS_TRANSFORMATIONS_DMN_2_TO_Y_LM_DMN_H

namespace TRAFOR
{
  /*!
   *  \class   TRANSFORM
   *  \ingroup TRANSFORM
   *
   *  \author  Peter Staar
   *  \brief   ...
   */
  template<typename type_input_0, typename type_input_1>
  struct TRANSFORM<dmn_2<dmn_0<type_input_0>, dmn_0<type_input_1> >, dmn_0<Y_lm_domain> >
  {
    const static bool VERBOSE = false;

    typedef Y_lm_domain output_type;
    
    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
			FUNC_LIB::function<scalartype_output, domain_output>& f_output,
			bool do_all_domains=false);
    
    template<typename scalartype, class domain_input, class domain_output>
    static void execute(FUNC_LIB::function<scalartype, domain_input> & f_input,
			FUNC_LIB::function<scalartype, domain_output>& f_output)
    {    
      typedef typename domain_input ::this_type type_list_input;
      typedef typename domain_output::this_type type_list_output;
    }

  private:

  };

}

#endif
