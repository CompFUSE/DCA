//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_DMN_0_TO_DMN_0_H
#define BASIS_TRANSFORMATIONS_DMN_0_TO_DMN_0_H

namespace MATH_ALGORITHMS
{
  /*!
   *  \class   TRANSFORM
   *  \ingroup TRANSFORM
   *
   *  \author  Peter Staar
   *  \brief   ...
   */
  template<typename type_input, typename type_output>
  class TRANSFORM<dmn_0<type_input>, dmn_0<type_output> >
  {
    const static bool VERBOSE = false;
    
  public:

    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
			FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {
      // typedef typename domain_input ::this_type type_list_input;
      // typedef typename domain_output::this_type type_list_output;
      
      if(VERBOSE)
	print_types(f_input, f_output);

      typedef typename SWAP_FIRST<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;      
      GENERIC_ASSERT<IS_EQUAL<TRANSFORMED_DOMAIN, domain_output>::CHECK>::execute();
      
      TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_first(f_input, f_output);
    }

    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output,
	     typename scalartype_T>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
			FUNC_LIB::function<scalartype_output, domain_output>& f_output,
			LIN_ALG::matrix<scalartype_T, LIN_ALG::CPU>& T)
    {
      // typedef typename domain_input ::this_type type_list_input;
      // typedef typename domain_output::this_type type_list_output;
      
      if(VERBOSE)
	print_types(f_input, f_output);

      typedef typename SWAP_FIRST<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;      
      GENERIC_ASSERT<IS_EQUAL<TRANSFORMED_DOMAIN, domain_output>::CHECK>::execute();
      
      TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_first(f_input, f_output, T);
    }

    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute_on_all(FUNC_LIB::function<scalartype_input, domain_input> & f_input,
			       FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {    
      // typedef typename domain_input ::this_type type_list_input;
      // typedef typename domain_output::this_type type_list_output;
      
      if(VERBOSE)
	print_types(f_input, f_output);

      typedef typename SWAP_ALL<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;

      GENERIC_ASSERT<IS_EQUAL<TRANSFORMED_DOMAIN, domain_output>::CHECK>::execute();
      
      TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_all(f_input, f_output);
    }

    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output,
	     typename scalartype_T>
    static void execute_on_all(FUNC_LIB::function<scalartype_input, domain_input> & f_input,
			       FUNC_LIB::function<scalartype_output, domain_output>& f_output,
			       LIN_ALG::matrix<scalartype_T, LIN_ALG::CPU>& T)
    {    
      // typedef typename domain_input ::this_type type_list_input;
      // typedef typename domain_output::this_type type_list_output;
      
      if(VERBOSE)
	print_types(f_input, f_output);

      typedef typename SWAP_ALL<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;

      GENERIC_ASSERT<IS_EQUAL<TRANSFORMED_DOMAIN, domain_output>::CHECK>::execute();
      
      TRANSFORM_DOMAINWISE<domain_input, domain_output, type_input, type_output>::execute_on_all(f_input, f_output, T);
    }

  private:
    
    template<typename scalartype_input , class domain_input, 
	     typename scalartype_output, class domain_output>
    static void print_types(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
			    FUNC_LIB::function<scalartype_output, domain_output>& f_output,
			    bool do_all_domains=false)
    {
      typedef typename domain_input ::this_type type_list_input;
      typedef typename domain_output::this_type type_list_output;
      
      print_type<type_input >::to_JSON(std::cout);
      std::cout << "\n\n";
      print_type<type_output>::to_JSON(std::cout);
      std::cout << "\n\n";
      
      if(do_all_domains)
	{
	  printTL<type_list_input >::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  typedef typename SWAP_ALL<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;
	  
	  printTL<typename TRANSFORMED_DOMAIN::this_type >::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  printTL<type_list_output>::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  FUNC_LIB::function<scalartype_output, TRANSFORMED_DOMAIN> T;
	  T       .print_fingerprint();
	  f_output.print_fingerprint();
	}
      else
	{
	  printTL<type_list_input >::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  typedef typename SWAP_FIRST<domain_input, type_input, type_output>::Result TRANSFORMED_DOMAIN;
	  
	  printTL<typename TRANSFORMED_DOMAIN::this_type >::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  printTL<type_list_output>::to_JSON(std::cout);
	  std::cout << "\n\n";
	  
	  FUNC_LIB::function<scalartype_output, TRANSFORMED_DOMAIN> T("FUNC_LIB::function<scalartype_output, TRANSFORMED_DOMAIN>");

	  f_input .print_fingerprint();
	  T       .print_fingerprint();
	  f_output.print_fingerprint();
	}
    }

  };

}

#endif
