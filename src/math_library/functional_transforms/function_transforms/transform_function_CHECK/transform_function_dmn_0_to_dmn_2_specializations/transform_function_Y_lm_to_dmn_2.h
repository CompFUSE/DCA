//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_Y_LM_DMN_TO_DMN_2_H
#define BASIS_TRANSFORMATIONS_Y_LM_DMN_TO_DMN_2_H

namespace TRAFOR
{
  /*!
   *  \class   TRANSFORM
   *  \ingroup TRANSFORM
   *
   *  \author  Peter Staar
   *  \brief   ...
   */
  template<typename type_output_0, typename type_output_1>
  struct TRANSFORM<dmn_0<Y_lm_domain>, dmn_2<dmn_0<type_output_0>, dmn_0<type_output_1> > >
  {
  public:

    typedef Y_lm_domain type_input;

    const static bool VERBOSE = false;

  public:
    
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
      
      if(VERBOSE)
	cout << __PRETTY_FUNCTION__ << endl;
      
      //typedef typename SWAP<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;
      
      
      //TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_first(f_input, f_output);
    }

  };
  
}

#endif
