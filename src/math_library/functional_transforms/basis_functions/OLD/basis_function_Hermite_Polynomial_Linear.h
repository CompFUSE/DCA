//-*-C++-*-

#ifndef EVALUATE_BASIS_FUNCTION_HERMITE_POLYNOMIAL_LINEAR_H
#define EVALUATE_BASIS_FUNCTION_HERMITE_POLYNOMIAL_LINEAR_H

namespace TRAFOR
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  class evaluate_basis_function<lh_dmn_type, rh_dmn_type, HERMITE_POLYNOMIAL_LINEAR>
  {
  public:

    const static bool ORTHOGONAL = false;

    const static BOUNDARY_CONDITIONS_type BC = lh_dmn_type::parameter_type::BOUNDARY_CONDITION;

    typedef typename lh_dmn_type::parameter_type::scalar_type scalartype;

  public:

    template<typename f_scalartype>
    static void execute(function<f_scalartype, lh_dmn_type>& f, int j)
    {
      scalartype y = rh_dmn_type::get_elements()[j];
      scalartype d = rh_dmn_type::get_elements()[1]-rh_dmn_type::get_elements()[0];

      for(int i=0; i<f.size(); i++)
	{
	  scalartype x = lh_dmn_type::get_elements()[i];

	  f(i) = spline(x, y, d);
	}
    }

    template<typename f_scalartype>
    static void norm(function<f_scalartype, rh_dmn_type>& f)
    {
      for(int i=0; i<f.size(); i++)
	f(i) = 1.;
    }

    template<typename f_scalartype>
    static void inner_product_function(function<f_scalartype, lh_dmn_type>& f)
    {

    }

    template<typename f_scalartype>
    static scalartype contraction(int i, int j)
    {
      
    }

  private:
    
    static scalartype spline(scalartype x, scalartype y, scalartype d)
    {
      scalartype delta = std::abs((y-x)/d);

      if(delta>=0-1.e-6 and delta<1)
	return 1.-delta;

      return 0;
    }
    
  };

}

#endif
