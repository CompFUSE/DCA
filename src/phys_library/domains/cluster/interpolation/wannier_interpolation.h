//-*-C++-*-

#ifndef WANNIER_INTERPOLATION_H
#define WANNIER_INTERPOLATION_H

/*! 
 *  \defgroup INTERPOLATION
 *  \ingroup  ALGORITHMS
 */

/*!
 *  \class   wannier_interpolation_domain_type
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   ...
 */
template<typename dmn_type, typename type_input, typename type_output>
struct wannier_interpolation_domain_type
{
//   typedef typename TL::Swap<dmn_type,type_input,type_output>::Result Result;

  typedef typename dmn_type::this_type dmn_type_list;
  typedef typename TL::Swap<dmn_type_list,type_input,type_output>::Result Result;
};

/*! \class   wannier_interpolation_kernel
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Wannier interpolation technique.
 */
template<typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation_kernel
{};

/*!
 *  \ingroup INTERPOLATION-TRANSFORM
 *
 *  \author  Peter Staar
 *  \brief   Perform a wannier-interpolation on a function.
 */
template<typename type_input, typename type_output, int dmn_number>
struct wannier_interpolation_any_2_any
{
  template<typename scalartype_1, typename dmn_type_1, typename scalartype_2, typename dmn_type_2>
  static void execute(FUNC_LIB::function<scalartype_1, dmn_type_1>& f_source, 
		      FUNC_LIB::function<scalartype_2, dmn_type_2>& f_target)
  {
    int Nb_sbdms    = f_source.signature();
    int Nb_elements = f_source.size();

    int* coordinate = new int[Nb_sbdms];
    memset(coordinate,0,sizeof(int)*Nb_sbdms);
    
    std::complex<double>* input_values  = new std::complex<double>[f_source[dmn_number] ];
    std::complex<double>* output_values = new std::complex<double>[f_target[dmn_number] ];

    {
      wannier_interpolation_kernel<type_input, type_output> kernel;
      
      int Nb_WI = Nb_elements/f_source[dmn_number];
      
      for(int l=0; l<Nb_WI; l++)
	{
	  int linind = l;
	  for(int j=Nb_sbdms-1; j>-1; j--)
	    {
	      if(j != dmn_number)
		{
		  coordinate[j] = linind % f_source[j];
		  linind = (linind-coordinate[j])/f_source[j];
		}
	    }
	  
	  f_source.slice(dmn_number, coordinate, input_values);

	  kernel.execute(input_values, output_values);

	  f_target.distribute(dmn_number, coordinate, output_values);
	}
    }

    delete [] coordinate;
    delete [] input_values;
    delete [] output_values;    
  }

};

/*!
 *  \class   wannier_interpolation_generic
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements the generic loop over all the subdomains.
 */
template<typename type_list1, typename type_list2, 
	 typename type_input, typename type_output, 
	 int dmn_shift, int next_index>
struct wannier_interpolation_generic
{
  template<typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input>& f_input, 
		      FUNC_LIB::function<scalartype_output, domain_output>& f_output)
  {
    //typedef typename TypeListAt<type_list1,IndexOf<type_list1, type_input>::value>::Result new_typelist1;
    //typedef typename TypeListAt<type_list2,IndexOf<type_list1, type_input>::value>::Result new_typelist2;

    wannier_interpolation_any_2_any<type_input, type_output, IndexOf<type_list1, type_input>::value + dmn_shift>::execute(f_input,f_output);

  }
};

template<typename type_list1, typename type_list2, typename type_input, typename type_output, int dmn_shift>
struct wannier_interpolation_generic<type_list1, type_list2, type_input, type_output, dmn_shift, -1>
{
  template<typename scalartype_1, typename dmn_type_1, typename scalartype_2, typename dmn_type_2>
  static void execute(FUNC_LIB::function<scalartype_1, dmn_type_1>& f_source, FUNC_LIB::function<scalartype_2, dmn_type_2>& F_target)
  {
    cout << "STOP" << endl;
  }
};

/*! \class   wannier_interpolation
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Wannier interpolation technique.
 */
template<typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation
{
  //#include "type_definitions.h"

public:
  
  template<typename scalartype_input, class domain_input, 
	   typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
		      FUNC_LIB::function<scalartype_output, domain_output>& f_output);
};

template<typename source_dmn_type, typename target_dmn_type>
template<typename scalartype_input, class domain_input, 
	 typename scalartype_output, class domain_output>
void wannier_interpolation<source_dmn_type, target_dmn_type>::execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
								      FUNC_LIB::function<scalartype_output, domain_output>& f_output)
{    
  typedef typename wannier_interpolation_domain_type<domain_input, source_dmn_type, target_dmn_type>::Result wannier_interpolation_domain;

//   GENERIC_ASSERT< IS_EQUAL_TYPE<domain_output, wannier_interpolation_domain>::check >::execute();

  typedef typename domain_output::this_type domain_output_list_type;
  GENERIC_ASSERT< IS_EQUAL_TYPE<domain_output_list_type, wannier_interpolation_domain>::check >::execute();
  
  typedef typename domain_input::this_type type_list_input;
  typedef typename domain_output::this_type type_list_output;
  
  GENERIC_ASSERT< (IndexOf<type_list_input, source_dmn_type>::value > -1) >::execute();
  
  wannier_interpolation_generic<type_list_input, 
    type_list_output, 
    source_dmn_type, 
    target_dmn_type, 
    0, 
    IndexOf<type_list_input, source_dmn_type>::value >::execute(f_input, f_output);
}


template<typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation<dmn_0<source_dmn_type>, dmn_0<target_dmn_type> >
{
public:

  template<typename scalartype_input, class domain_input, 
	   typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
		      FUNC_LIB::function<scalartype_output, domain_output>& f_output)
  {
    wannier_interpolation<source_dmn_type, target_dmn_type>::execute(f_input,f_output);
  }
};

#endif
