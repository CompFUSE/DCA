//-*-C++-*-

#ifndef COVARIANCE_FUNCTION_H
#define COVARIANCE_FUNCTION_H

namespace MATH_LIBRARY
{
  enum GP_COVARIANCE {LINEAR_EXPONENTIAL, PERIODIC_LINEAR_EXPONENTIAL,
                      SQUARED_EXPONENTIAL, PERIODIC_SQUARED_EXPONENTIAL};

  /*!
   *
   *    page 19, eqn 2.31, Rasmussen and Williams
   */
  template<GP_COVARIANCE COVARIANCE, typename k_dmn_t>
  class covariance_function
  {
    const static int DIMENSION = k_dmn_t::parameter_type::DIMENSION;

  public:

    covariance_function();
    ~covariance_function();

    double execute(std::vector<double>& x_i);
  };

}

#endif
