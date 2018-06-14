// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
