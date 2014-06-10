//-*-C++-*-

#ifndef BASIS_FUNCTION_HARMONICS_FERMIONIC_POSITIVE_H
#define BASIS_FUNCTION_HARMONICS_FERMIONIC_POSITIVE_H

namespace TRAFOR
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  class basis_function<lh_dmn_type, rh_dmn_type, HARMONICS_FERMIONIC_POSITIVE>
  {
  public:

    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;
    
    typedef std::complex<lh_scalar_type> f_scalar_type;

  public:

    static f_scalar_type execute(int i, int j)
    {
      return execute(lh_dmn_type::get_elements()[i], rh_dmn_type::get_elements()[j]);
    }

    static f_scalar_type execute(lh_element_type& lh_elem, rh_element_type& rh_elem)
    {
      const static f_scalar_type I(0,1);

      f_scalar_type phase = dot_prod(lh_elem, rh_elem);

      return std::exp(I*phase);
    }

  private:
    
    template<typename scalartype>
    inline static scalartype dot_prod(scalartype x, scalartype y)
    {
      return x*y;
    }

    template<typename scalartype>
    inline static scalartype dot_prod(std::vector<scalartype> x, std::vector<scalartype> y)
    {
      scalartype result=0;
      for(size_t l=0; l<x.size(); l++)
	result += x[l]*y[l];
      return result;
    }
  };

}

#endif
