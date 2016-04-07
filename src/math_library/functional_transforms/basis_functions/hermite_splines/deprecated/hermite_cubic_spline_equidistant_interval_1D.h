//-*-C++-*-

#ifndef BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_INTERVAL_1D_H
#define BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_INTERVAL_1D_H

namespace math_algorithms
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  struct hermite_cubic_spline<lh_dmn_type, rh_dmn_type, EQUIDISTANT, INTERVAL, 1>
  {
    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;
    
    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;
    
    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;
    
    typedef lh_scalar_type f_scalar_type;
    
    static f_scalar_type execute(int i, int j)
    {
      lh_scalar_type x = lh_dmn_type::get_elements()[i];
      rh_scalar_type y = rh_dmn_type::get_elements()[j];
      
      const static rh_scalar_type a = -0.5;
      lh_scalar_type              d = rh_dmn_type::get_elements()[1]-rh_dmn_type::get_elements()[0];
      
      rh_scalar_type delta = std::abs((y-x)/d);
      
      f_scalar_type result = 0;
      
      if(delta<2)
	{
	  result = hermite_spline::cubic(x, y, d, a);
	  
	  int N               = rh_dmn_type::get_size();
	  rh_element_type min = rh_dmn_type::get_elements()[0];
	  rh_element_type max = rh_dmn_type::get_elements()[N-1];
	  
	  if(j==0)
	    result += hermite_spline::cubic(x, min-d, d, a)*(3.);
	  
	  if(j==1)
	    result -= hermite_spline::cubic(x, min-d, d, a)*(3.);
	  
	  if(j==2)
	    result += hermite_spline::cubic(x, min-d, d, a)*(1.);
	  
	  if(j==N-1)
	      result += hermite_spline::cubic(x, max+d, d, a)*(3.);
	  
	  if(j==N-2)
	    result -= hermite_spline::cubic(x, max+d, d, a)*(3.);
	  
	  if(j==N-3)
	    result += hermite_spline::cubic(x, max+d, d, a)*(1.);
	}
      
      return result;
    }
    
  };
  
}

#endif
