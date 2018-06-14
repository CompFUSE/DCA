// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements numerical_error_domain.hpp.

#include "dca/phys/domains/quantum/numerical_error_domain.hpp"
#include <cmath>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::vector<double>& numerical_error_domain::initialize_elements() {
  static std::vector<element_type> v(0);

  for (int i = -16; i < 0; i++)
    for (double j = 1; j < 10; j += 1.)
      v.push_back(j * std::pow(10., i));

  return v;
}

}  // domains
}  // phys
}  // dca
