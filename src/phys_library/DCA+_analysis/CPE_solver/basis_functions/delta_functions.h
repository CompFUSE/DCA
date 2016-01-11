//-*-C++-*-

#ifndef DELTA_FUNCTION_H
#define DELTA_FUNCTION_H

/*!
 *  \ingroup CPE-BASIS-FUNCTIONS
 *
 *  \author  Peter Staar
 */
class delta_function
{
#include "type_definitions.h"

public:

  typedef double element_type;

public:

  static int&                  get_size();
  static std::vector<double>& get_elements();

  template<typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double& epsilon();

  static std::complex<double> phi(int n, std::complex<double> z);
};

template<typename parameters_type>
void delta_function::initialize(parameters_type& parameters)
{
  cout << __PRETTY_FUNCTION__ << endl;

  get_size()     = parameters.get_poles().size();
  get_elements() = parameters.get_poles();

  epsilon()      = 2.*parameters.get_EPSILON();
}

int& delta_function::get_size()
{
  static int size = 0;
  return size;
}

std::vector<double>& delta_function::get_elements()
{
  static std::vector<double> elements(0,0);
  return elements;
}

double& delta_function::epsilon()
{
  static double epsilon = 1.e-3;
  return epsilon;
}

std::complex<double> delta_function::phi(int n, std::complex<double> z)
{
  assert(n>=0 && n<get_size());

  double x = get_elements()[n];

  std::complex<double> A_mn;

  if(imag(z) > epsilon()){
    real(A_mn) = real(1./(z-x));
    imag(A_mn) = imag(1./(z-x));
  }
  else{
    real(A_mn) = std::fabs(real(z)-x)<1.e-6?  0. : 1/(real(z)-x);
    imag(A_mn) = 0.;
  }

  return A_mn;
}

#endif
