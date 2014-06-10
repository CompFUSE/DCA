//-*-C++-*-

#ifndef BASIS_FUNCTION_LEGENDRE_LM_H
#define BASIS_FUNCTION_LEGENDRE_LM_H

namespace TRAFOR
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  class basis_function<lh_dmn_type, rh_dmn_type, LEGENDRE_LM>
  {
  public:

    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;
    
    typedef lh_scalar_type f_scalar_type;

  public:

    static f_scalar_type execute(int i, int j)
    {
      return execute(lh_dmn_type::get_elements()[i], rh_dmn_type::get_elements()[j]);
    }

    static f_scalar_type execute(lh_element_type& lh_elem, rh_element_type& rh_elem)
    {
      lh_scalar_type delta = get_delta();
      lh_scalar_type min   = lh_dmn_type::get_elements()[0];

      return compute(2*(lh_elem-min)/delta-1, rh_elem);
    }

  private:

    inline static lh_scalar_type get_delta()
    {
      lh_scalar_type result = 0;

      switch(lh_spec_dmn_type::boundary_condition)
 	{
 	case INTERVAL:
 	  result = (lh_dmn_type::get_elements()[1]-lh_dmn_type::get_elements()[0])*(lh_dmn_type::get_size()-1);
 	  break;

	case PERIODIC:
	case ANTIPERIODIC:
	  result = (lh_dmn_type::get_elements()[1]-lh_dmn_type::get_elements()[0])*(lh_dmn_type::get_size()+0);
	  break;

	default:
	  throw std::logic_error(__FUNCTION__);
	}

      return result;
    }

    inline static f_scalar_type compute(lh_scalar_type t, std::pair<int,int> i)
    {
      assert(t>-1.-1.e-6 and t<1.+1.e-6);

      int l = i.first;
      int m = i.second;

      f_scalar_type result = compute_positive(t, l, abs_value(m));

      if(m<0)
	result *= std::pow(-1., m)*factorial(l-m)/factorial(l+m);
	  
      return result;
    }

    inline static f_scalar_type compute_positive(lh_scalar_type t, int l, int m)
    {
      assert(t>-1.-1.e-6 and t<1.+1.e-6);

      f_scalar_type result = 0;

      if(m>l)
	result = 0.;

      if(l==0 and m==0)
	return 1.;

      if(l==m and t>-1.+1.e-6 and t<1.-1.e-6)
	result = std::pow(-1., l)*double_factorial(2*l-1)*std::pow(1.-t*t, l/2.);
      else
	result = 0;

      if(m<l)
	result = 1./f_scalar_type(l-m)*((2*l-1)*t*compute_positive(t,l-1,m)-(l-1+m)*compute_positive(t,l-2,m));

      assert(result==result);

      return result;
    }    

    inline static int abs_value(int m)
    {
      return m<0? -m : m;
    }

    inline static int factorial(int m)
    {
      assert(m>-1);

      switch(m)
	{
	case 0:
	case 1:
	  return 1; 

	default:
	  return m*factorial(m-1);
	}
    }

    inline static int double_factorial(int m)
    {
      assert(m>-1);
      assert((m+1)%2==0);

      switch(m)
	{
	case 0:
	case 1:
	  return 1; 

	default:
	  return m*factorial(m-2);
	}
    }

  };

}

#endif
