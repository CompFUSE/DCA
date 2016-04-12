//-*-C++-*-

#ifndef HSPLINE_INTERPOLATION_ALGORITHM_H
#define HSPLINE_INTERPOLATION_ALGORITHM_H

#include <include/dca/util/type_utils.hpp>

/*!
 *  \class   hspline_interpolation_domain_type
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   ...
 */
template<typename dmn_type, typename type_input, typename type_output>
struct hspline_interpolation_domain_type
{
  typedef typename dmn_type::this_type dmn_type_list;
  typedef typename dca::util::Swap<dmn_type_list,type_input,type_output>::type Result;
};

/*! \class   hspline_interpolation_kernel
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Hspline interpolation technique.
 */
template<typename scalartype, typename source_dmn_type, typename target_dmn_type>
class hspline_interpolation_kernel
{};

/*!
 *  \ingroup INTERPOLATION-TRANSFORM
 *
 *  \author  Peter Staar
 *  \brief   Perform a hspline-interpolation on a function.
 */
template<typename type_input, typename type_output, int dmn_number>
struct hspline_interpolation_any_2_any
{
  template<typename scalartype, typename dmn_type_1, typename dmn_type_2>
  static void execute(FUNC_LIB::function<scalartype, dmn_type_1>& f_source,
                      FUNC_LIB::function<scalartype, dmn_type_2>& f_target,
                      double a)
  {
    int Nb_sbdms    = f_source.signature();
    int Nb_elements = f_source.size();

    int* coordinate = new int[Nb_sbdms];
    memset(coordinate,0,sizeof(int)*Nb_sbdms);

    scalartype* input_values  = new scalartype[f_source[dmn_number] ];
    scalartype* output_values = new scalartype[f_target[dmn_number] ];

    {
      hspline_interpolation_kernel<scalartype, type_input, type_output> kernel(a);

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
 *  \class   hspline_interpolation_generic
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements the generic loop over all the subdomains.
 */
template<typename type_list1, typename type_list2,
         typename type_input, typename type_output,
         int dmn_shift, int next_index>
struct hspline_interpolation_generic
{
  template<typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input>& f_input,
                      FUNC_LIB::function<scalartype_output, domain_output>& f_output,
                      double a)
  {
    //typedef typename TypeListAt<type_list1,IndexOf<type_list1, type_input>::value>::Result new_typelist1;
    //typedef typename TypeListAt<type_list2,IndexOf<type_list1, type_input>::value>::Result new_typelist2;

    hspline_interpolation_any_2_any<type_input, type_output, IndexOf<type_input, type_list1>::value + dmn_shift>::execute(f_input,f_output, a);
  }
};

template<typename type_list1, typename type_list2, typename type_input, typename type_output, int dmn_shift>
struct hspline_interpolation_generic<type_list1, type_list2, type_input, type_output, dmn_shift, -1>
{
  template<typename scalartype_1, typename dmn_type_1, typename scalartype_2, typename dmn_type_2>
    static void execute(FUNC_LIB::function<scalartype_1, dmn_type_1>& /*f_source*/, FUNC_LIB::function<scalartype_2, dmn_type_2>& /*F_target*/, double /*a*/)
  {}
};




/*! \class   hspline_interpolation
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Hspline interpolation technique.
 */
template<typename source_dmn_type, typename target_dmn_type>
class hspline_interpolation
{

public:

  template<typename scalartype_input, class domain_input,
           typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
                      FUNC_LIB::function<scalartype_output, domain_output>& f_output,
                      double a)
  {

    typedef typename hspline_interpolation_domain_type<domain_input, source_dmn_type, target_dmn_type>::Result hspline_interpolation_domain;

    typedef typename domain_output::this_type domain_output_list_type;
    dca::util::assert_same<domain_output_list_type, hspline_interpolation_domain>();

    typedef typename domain_input::this_type type_list_input;
    typedef typename domain_output::this_type type_list_output;

    static_assert(IndexOf<source_dmn_type, type_list_input>::value > -1, "Type list error");

    hspline_interpolation_generic<type_list_input,
                                  type_list_output,
                                  source_dmn_type,
                                  target_dmn_type,
                                  0,
                                  IndexOf<source_dmn_type, type_list_input>::value >::execute(f_input, f_output, a);
  }

  /*
    template<typename scalartype_input, class domain_input,
    typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<std::complex<scalartype_input>, domain_input> & f_input,
    FUNC_LIB::function<std::complex<scalartype_output>, domain_output>& f_output,
    double a);
  */
};











template<typename source_dmn_type, typename target_dmn_type>
class hspline_interpolation<dmn_0<source_dmn_type>, dmn_0<target_dmn_type> >
{
public:

  template<typename scalartype_input, class domain_input,
           typename scalartype_output, class domain_output>
  static void execute(FUNC_LIB::function<scalartype_input , domain_input> & f_input,
                      FUNC_LIB::function<scalartype_output, domain_output>& f_output,
                      double a)
  {
    hspline_interpolation<source_dmn_type, target_dmn_type>::execute(f_input,f_output,a);
  }
};

#endif
