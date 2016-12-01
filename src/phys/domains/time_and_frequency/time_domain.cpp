// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements time_domain.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include <cmath>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

void time_domain::initialize(const double beta, const int time_slices, const double eps) {
  get_size() = 2 * (time_slices + 1);
  get_elements().resize(get_size());

  for (int i = 0; i < get_size() / 2; i++) {
    get_elements()[i + get_size() / 2] = double(i) / (double(get_size()) / 2. - 1.) * beta;
    get_elements()[i] = -beta + double(i) / (double(get_size()) / 2. - 1.) * beta;
  }

  get_elements()[0] += eps;
  get_elements()[get_size() - 1] -= eps;
  get_elements()[get_size() / 2] += eps;
  get_elements()[get_size() / 2 - 1] -= eps;
}

void time_domain::initialize_integration_domain(int level, std::vector<scalar_type>& weights,
                                                std::vector<element_type>& elements) {
  int N = std::pow(2., level);

  weights.resize(0);
  elements.resize(0);

  for (int i = 0; i < get_size() / 2 - 1; i++) {
    element_type min = get_elements()[i];
    element_type max = get_elements()[i + 1];

    double delta = double(max - min) / double(N);

    for (int j = 0; j < N; j++) {
      elements.push_back(min + j * delta);

      weights.push_back(get_volume());
    }
  }

  for (int i = get_size() / 2; i < get_size() - 1; i++) {
    element_type min = get_elements()[i];
    element_type max = get_elements()[i + 1];

    double delta = double(max - min) / double(N);

    for (int j = 0; j < N; j++) {
      elements.push_back(min + j * delta);

      weights.push_back(get_volume());
    }
  }
}

}  // domains
}  // phys
}  // dca
