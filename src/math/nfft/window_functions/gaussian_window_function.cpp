// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements gaussian_window_function.hpp

#include "dca/math/nfft/window_functions/gaussian_window_function.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft

int gaussian_window_function::n = 1;
int gaussian_window_function::m = 1;
double gaussian_window_function::sigma = 1;

}  // nfft
}  // math
}  // dca
