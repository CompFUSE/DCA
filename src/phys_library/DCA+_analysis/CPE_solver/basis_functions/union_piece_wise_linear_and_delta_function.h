//-*-C++-*-

#ifndef UNION_PIECE_WISE_LINEAR_AND_DELTA_FUNCTION_H
#define UNION_PIECE_WISE_LINEAR_AND_DELTA_FUNCTION_H

/*!
 *  \ingroup CPE-BASIS-FUNCTIONS
 *
 *  \author  Peter Staar
 */
class union_piece_wise_linear_and_delta_function
{
#include "type_definitions.h"

public:

  typedef double element_type;

public:

  static int&                 get_size();
  static std::vector<double>& get_elements();

  template<typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n); 

  static std::complex<double> phi(int n, std::complex<double> z);
};

int& union_piece_wise_linear_and_delta_function::get_size()
{
  static int size = 0;
  return size;
}

std::vector<double>& union_piece_wise_linear_and_delta_function::get_elements()
{
  static std::vector<double> elements(0);
  return elements;
}

template<typename parameters_type>
void union_piece_wise_linear_and_delta_function::initialize(parameters_type& parameters)
{
  delta_function            ::initialize(parameters);
  piece_wise_linear_function::initialize(parameters);

  get_size() = piece_wise_linear_function::get_size() + delta_function::get_size();

  get_elements().insert(get_elements().end(), piece_wise_linear_function::get_elements().begin(), piece_wise_linear_function::get_elements().end());
  get_elements().insert(get_elements().end(), delta_function::get_elements().begin()            , delta_function::get_elements().end());

  cout << __PRETTY_FUNCTION__ 
       << "\n\t" <<  piece_wise_linear_function::get_size() 
       << "\t" <<  delta_function::get_size() 
       << "\t" <<  get_size() 
       << "\n";
}

double union_piece_wise_linear_and_delta_function::volume(int n)
{
  assert(n>=0 && n<get_size());

  return piece_wise_linear_function::volume(0);
}

std::complex<double> union_piece_wise_linear_and_delta_function::phi(int n, std::complex<double> z)
{
  assert(n>=0 && n<get_size());

  if(n<piece_wise_linear_function::get_size())
    return piece_wise_linear_function::phi(n,z);
  else
    return delta_function::phi(n-piece_wise_linear_function::get_size(),z);
}

#endif
