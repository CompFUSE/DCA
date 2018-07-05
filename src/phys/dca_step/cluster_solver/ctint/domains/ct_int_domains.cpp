// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include <algorithm>

#include "dca/phys/dca_step/cluster_solver/ctint/domains/ct_int_domains.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {

void PositiveTimeDomain::initialize() {
  if (not TimeDomain::is_initialized())
    throw(std::logic_error("Time domain was not initialized."));
  if (initialized_) {  // Nothing to do (assume time domain did not change)
    assert(elements_.size() == TimeDomain::get_size() / 2);
    return;
  }

  const std::size_t size = TimeDomain::get_size() / 2;
  elements_.resize(size);
  std::copy_n(&TimeDomain::get_elements()[size], size, elements_.begin());

  // TODO consider if  you have to remove the epsilon shift
  // Remove shift by epsilon
  elements_[0] = 0.;
  elements_.back() = TimeDomain::get_volume();  // = beta.
  initialized_ = true;
}

bool PositiveTimeDomain::initialized_ = false;
std::string PositiveTimeDomain::name_ = "positive time domain";
std::vector<typename PositiveTimeDomain::ElementType> PositiveTimeDomain::elements_;

int G0ParametersDomain::size_ = -1;

}  // dca
}  // phys
}  // solver
}  // ctint
