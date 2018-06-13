// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements time_domain.hpp.

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include <cmath>  // for std::pow
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

bool time_domain::initialized_ = false;
const std::string time_domain::name_ = "time-domain";
time_domain::scalar_type time_domain::beta_ = -1.;
std::vector<time_domain::element_type> time_domain::elements_;

void time_domain::initialize(const scalar_type beta, const int time_slices, const scalar_type eps) {
  if (initialized_)
    throw std::logic_error("time_domain has already been initialized.");

  beta_ = beta;

  const int size = 2 * (time_slices + 1);
  elements_.resize(size);

  const element_type step = beta / time_slices;

  for (std::size_t i = 0; i < size / 2; i++) {
    elements_[i + size / 2] = i * step;
    elements_[i] = -beta + i * step;
  }

  // Shift boundaries and zeros by a small number 'eps'.
  elements_[0] += eps;
  elements_[size - 1] -= eps;
  elements_[size / 2] += eps;
  elements_[size / 2 - 1] -= eps;

  initialized_ = true;
}

void time_domain::initialize_integration_domain(const int level, std::vector<scalar_type>& weights,
                                                std::vector<element_type>& nodes) {
  if (!initialized_)
    throw std::logic_error("time_domain has not been initialized.");

  const int num_nodes_per_time_slice = std::pow(2., level);

  weights.clear();
  nodes.clear();

  for (std::size_t i = 0; i < elements_.size() / 2 - 1; i++) {
    const element_type min = elements_[i];
    const element_type max = elements_[i + 1];
    const element_type delta = (max - min) / num_nodes_per_time_slice;

    for (int j = 0; j < num_nodes_per_time_slice; j++) {
      nodes.push_back(min + j * delta);
      weights.push_back(get_volume());
    }
  }

  for (std::size_t i = elements_.size() / 2; i < elements_.size() - 1; i++) {
    const element_type min = elements_[i];
    const element_type max = elements_[i + 1];
    const element_type delta = (max - min) / num_nodes_per_time_slice;

    for (int j = 0; j < num_nodes_per_time_slice; j++) {
      nodes.push_back(min + j * delta);
      weights.push_back(get_volume());
    }
  }
}

}  // domains
}  // phys
}  // dca
