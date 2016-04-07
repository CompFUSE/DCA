//-*-C++-*-

#ifndef UNION_RADIAL_AND_DELTA_FUNCTION_H
#define UNION_RADIAL_AND_DELTA_FUNCTION_H
#include"phys_library/domain_types.hpp"
using namespace types;

/*!
 *  \ingroup CPE-BASIS-FUNCTIONS
 *
 *  \author  Peter Staar
 */
class union_radial_and_delta_function
{

public:

  typedef double element_type;

public:

  static int                  get_size();
  static std::vector<double>& get_elements();

  template<typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double               volume(int n);
  static std::complex<double> phi(int n, std::complex<double> z);

private:
  
  static std::vector<double> poles;
};

std::vector<double> union_radial_and_delta_function::poles(0,0);

int union_radial_and_delta_function::get_size()
{
  return get_elements().size();
}

std::vector<double>& union_radial_and_delta_function::get_elements()
{
  static std::vector<double> elements(0);
  return elements;
}

template<typename parameters_type>
void union_radial_and_delta_function::initialize(parameters_type& parameters)
{
  poles = parameters.get_poles();
  
  get_elements() = w_REAL::get_elements();
}

double union_radial_and_delta_function::volume(int n)
{
  assert(n>=0 && n<w_REAL::dmn_size());

  double volume;

  for(size_t l1=0; l1<poles.size(); l1++)
    for(size_t l2=0; l2<get_elements().size(); l2++)
      if(std::fabs(poles[l1]-get_elements()[l2])<1.e-6)
	return M_PI;

  if(n==0){
    volume = 2.*(get_elements()[1]-get_elements()[0])/2.;
  }
  else{
    if(n==get_size()-1)
      volume = 2.*(get_elements()[n]-get_elements()[n-1])/2.;
    else
      volume = (get_elements()[n+1]-get_elements()[n-1])/2.;
  }

  return volume;
}

std::complex<double> union_radial_and_delta_function::phi(int n, std::complex<double> z)
{
  assert(n>=0 && n<get_size());

  double w=1.;

  std::complex<double> A_mn, x0, x1, x2;

  for(size_t l=0; l<poles.size(); l++){ // go over all poles
    if(abs(get_elements()[n]-poles[l])<1.e-6) // elements()[n] is a pole !!
      { 
	if(abs(imag(z)<1.e-6)) // z lies on the real axis
	  {
	    if(abs(z-poles[l])<1.e-6) // z is on the pole !
	      return 0.;
	    else
	      return std::complex<double>(w,0.)/(real(z)-poles[l]);
	  }
	else
	  return w/(z-poles[l]);
      }
  }

  if(n==0)
    {
      double delta_x = (get_elements()[1]-get_elements()[0]);

      x0 = get_elements()[0]-delta_x;
      x1 = get_elements()[0];
      x2 = get_elements()[0]+delta_x;
    }
  else
    {
      if(n==get_size()-1)
	{
	  double delta_x = (get_elements()[n]-get_elements()[n-1]);

	  x0 = get_elements()[n]-delta_x;
	  x1 = get_elements()[n];
	  x2 = get_elements()[n]+delta_x;
	}
      else
	{
	  x0 = get_elements()[n-1];
	  x1 = get_elements()[n];
	  x2 = get_elements()[n+1];
	}
    }

  real(A_mn) = real((-x0 + x1 + (x0 - z)*(std::log(x0 - z) - std::log(x1 - z)))/(x0 - x1) + (x1 - x2 - (x2 - z)*(std::log(-x1 + z) - std::log(-x2 + z)))/(x1 - x2));
  imag(A_mn) = imag((-x0 + x1 + (x0 - z)*(std::log(x0 - z) - std::log(x1 - z)))/(x0 - x1) + (x1 - x2 - (x2 - z)*(std::log(-x1 + z) - std::log(-x2 + z)))/(x1 - x2));

  assert(A_mn==A_mn); // no nan's !

  return A_mn;
}

#endif
