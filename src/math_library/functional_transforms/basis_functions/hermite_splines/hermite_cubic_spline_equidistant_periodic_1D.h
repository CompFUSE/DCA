//-*-C++-*-

#ifndef BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_PERIODIC_1D_H
#define BASIS_FUNCTION_HERMITE_CUBIC_SPLINE_EQUIDISTANT_PERIODIC_1D_H

namespace MATH_ALGORITHMS
{
  template<typename lh_dmn_type, typename rh_dmn_type>
  class hermite_cubic_spline<lh_dmn_type, rh_dmn_type, EQUIDISTANT, PERIODIC, 1>
  {
  private:

    typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
    typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

    typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
    typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

    typedef typename lh_spec_dmn_type::element_type lh_element_type;
    typedef typename rh_spec_dmn_type::element_type rh_element_type;

    typedef lh_scalar_type f_scalar_type;

  public:

    /*
      static f_scalar_type execute(int i, int j)
      {
      const static rh_scalar_type a = -0.5;

      lh_element_type x = lh_dmn_type::get_elements()[i];
      rh_element_type y = rh_dmn_type::get_elements()[j];

      int N            = rh_dmn_type::get_size();
      lh_scalar_type d = rh_dmn_type::get_elements()[1]-rh_dmn_type::get_elements()[0];
      rh_scalar_type D = d*N;

      f_scalar_type result = 0;
      for(int l0=-2; l0<=2; l0++){
      rh_scalar_type delta = std::abs((y-l0*D-x)/d);
      if(delta<2)
      result += hermite_spline::cubic(x, y-l0*D, d, a);
      }

      return result;
      }
    */

    static f_scalar_type execute(int i, int j)
    {
      const static rh_scalar_type a = -0.5;

      lh_element_type x = lh_dmn_type::get_elements()[i];
      rh_element_type y = rh_dmn_type::get_elements()[j];

      rh_scalar_type* basis       = rh_dmn_type::get_basis();
      rh_scalar_type* super_basis = rh_dmn_type::get_super_basis();

      rh_scalar_type* inv_basis       = rh_dmn_type::get_inverse_basis();
      rh_scalar_type* inv_super_basis = rh_dmn_type::get_inverse_super_basis();

      rh_element_type delta        = 0.;
      rh_element_type delta_affine = 0.;

      f_scalar_type result = 0;
      for(int l0=-2; l0<=2; l0++){

        delta = (y-x)-(l0*super_basis[0]);

        delta_affine += inv_super_basis[0]*delta;

        while(delta_affine>0.5-1.e-6)
          delta_affine -= 1.;

        while(delta_affine<-0.5-1.e-6)
          delta_affine += 1.;

        delta += super_basis[0]*delta_affine;

        delta_affine[0] += inv_basis[0]*delta;

        f_scalar_type t_result = (delta_affine>-2. and delta_affine<2.) ? hermite_spline::cubic(0., delta_affine, 1., a) : 0. ;

        result += t_result;
      }

      return result;
    }

  };

}

#endif
