// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This domain contained in this file describes the momentum exchange inside a two particle Green's
// Function.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_CLUSTER_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_CLUSTER_DOMAIN_HPP

#include <cassert>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class MomentumExchangeDomain {
public:
  using scalar_type = int;
  using element_type = int;

  // Initialize size and elements.
  template <class Parameters>
  static void initialize(const Parameters& parameters);

  // Returns 1 if Parameters::compute_all_transfers() == true, otherwise the cluster size
  static int get_size() {
    assert(initialized_);
    return elements_.size();
  }

  // Contains all indices relative to the exchanged momentum inside the two particle
  // computation.
  // Precondition: the domain is initialized.
  static inline const std::vector<int>& get_elements() {
    assert(initialized_);
    return elements_;
  }

  static const std::string& get_name() {
    const static std::string name = "Momentum exchange domain.";
    return name;
  }

  template <class Writer>
  static void write(Writer& writer);

private:
  static std::vector<int> elements_;
  static bool initialized_;
};
std::vector<int> MomentumExchangeDomain::elements_;
bool MomentumExchangeDomain::initialized_ = false;

template <class Parameters>
void MomentumExchangeDomain::initialize(const Parameters& parameters) {
  if (parameters.compute_all_transfers()) {
    const int size = Parameters::KClusterDmn::dmn_size();
    elements_.resize(size);
    int idx_value = 0;
    for (int& elem : elements_)
      elem = idx_value++;
  }

  else {
    elements_ = std::vector<int>{parameters.get_four_point_momentum_transfer_index()};
  }

  initialized_ = true;
}

template <class Writer>
void MomentumExchangeDomain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("element_indices_", elements_);
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_CLUSTER_DOMAIN_HPP
