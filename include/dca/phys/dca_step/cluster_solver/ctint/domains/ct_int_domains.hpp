// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides the POsitiveTimeDomains and ParametersDomain used by the g0 interpolator.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_POSITIVE_TIME_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_POSITIVE_TIME_DOMAIN_HPP

#include <stdexcept>

#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

class PositiveTimeDomain {
private:
  using TimeDomain = dca::phys::domains::time_domain;

public:
  using ElementType = typename TimeDomain::element_type;
  using element_type = ElementType;  // For putting PositiveTimeDomain in dmn_0.

  // Creates a domain that contains the elements in the range [0, beta] of the base domain.
  static void initialize();

  static std::size_t get_size() {
    return elements_.size();
  }
  static const std::string& get_name() {
    return name_;
  }
  static std::vector<ElementType>& get_elements() {
    assert(initialized_);
    return elements_;
  }

  static bool is_initialized() {
    return initialized_;
  }

private:
  static void commonInit(const std::string& name);

  static bool initialized_;
  static std::string name_;
  static std::vector<ElementType> elements_;
};

struct G0ParametersDomain {
  void static initialize(int n) {
    size_ = n;
  }

  void static reset() {
    size_ = -1;
  }

  // Used by dmn_variadic to compute the size.
  int static get_size() {
    return size_;
  }

  using element_type = void;
  using this_type = G0ParametersDomain;

private:
  static int size_;
};



}  // dca
}  // phys
}  // solver
}  // ctint

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_POSITIVE_TIME_DOMAIN_HPP
