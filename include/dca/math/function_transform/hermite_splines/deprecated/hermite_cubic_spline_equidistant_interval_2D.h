//-*-C++-*-

#ifndef BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_INTERVAL_2D_H
#define BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_INTERVAL_2D_H

namespace math_algorithms
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  struct hermite_cubic_spline<lh_dmn_type, rh_dmn_type, EQUIDISTANT, INTERVAL, 2>
  {
    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;
    
    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;
    
    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;
    
    typedef lh_scalar_type f_scalar_type;
    
    static f_scalar_type execute(int /*i*/, int /*j*/)
    {
      throw std::logic_error(__FUNCTION__);
    }
    
  };
  
}

#endif
