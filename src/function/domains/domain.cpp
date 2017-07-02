// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements domain.hpp.

#include "dca/function/domains/domain.hpp"

#include <cassert>
#include <cstdlib>

namespace dca {
namespace func {
// dca::func::

domain::domain()
    : size(0),
      branch_domain_sizes(0),
      leaf_domain_sizes(0),

      branch_domain_steps(0),
      leaf_domain_steps(0) {}

void domain::reset() {
  size = 0;

  leaf_domain_sizes.resize(0);
  branch_domain_sizes.resize(0);

  leaf_domain_steps.resize(0);
  branch_domain_steps.resize(0);
}

void domain::linind_2_subind(int linind, int* subind) const {
  assert(linind >= 0 && linind < size);

  for (std::size_t i = 0; i < leaf_domain_sizes.size(); i++) {
    subind[i] = linind % leaf_domain_sizes[i];
    linind = (linind - subind[i]) / leaf_domain_sizes[i];
  }

  assert(linind == 0);
}

void domain::subind_2_linind(const int* const subind, int& linind) const {
  linind = 0;

  for (int i = leaf_domain_sizes.size() - 1; i >= 0; i--) {
    assert(subind[i] < leaf_domain_sizes[i]);
    linind = subind[i] + linind * leaf_domain_sizes[i];
  }

  assert(linind >= 0 && linind < size);
}

}  // func
}  // dca
