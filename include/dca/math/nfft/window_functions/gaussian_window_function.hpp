// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a Gaussian window function.

#ifndef DCA_MATH_NFFT_WINDOW_FUNCTIONS_GAUSSIAN_WINDOW_FUNCTION_HPP
#define DCA_MATH_NFFT_WINDOW_FUNCTIONS_GAUSSIAN_WINDOW_FUNCTION_HPP

#include <cmath>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft

class gaussian_window_function {
public:
  static int n;
  static int m;
  static double sigma;

  static double phi_t(double x) {
    return 1. / std::sqrt(M_PI * b_val()) * std::exp(-(n * n * x * x) / b_val());
  }

  static double d_phi_t(double x) {
    return 1. / std::sqrt(M_PI * b_val()) * (-2. * n * n * x / b_val()) *
           std::exp(-(n * n * x * x) / b_val());
  }

  static double phi_wn(int wn) {
    return 1. / double(n) * std::exp(-b_val() * std::pow(M_PI * double(wn) / double(n), 2));
  }

private:
  static double b_val() {
    return m / M_PI * (2 * sigma) / (2 * sigma - 1);
  }
};

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_WINDOW_FUNCTIONS_GAUSSIAN_WINDOW_FUNCTION_HPP
