// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// The domain contained in this file describes the momentum exchange inside a two particle Green's
// Function.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_CLUSTER_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_CLUSTER_DOMAIN_HPP

#include <array>
#include <cassert>
#include <string>
#include <vector>

#include "dca/phys/domains/cluster/cluster_symmetry.hpp"

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
  static void initialize(const Parameters& parameters, bool verbose = false);

  // Returns the number of computed momentum exchanges.
  // Precondition: the domain is initialized.
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

  static bool isInitialized() {
    return initialized_;
  }

  template <class Writer>
  static void write(Writer& writer);

private:
  static void initialize(bool compute_all_transfers, int transfer_index,
                         const std::vector<std::vector<double>>& cluster_elements,
                         const std::vector<std::array<int, 2>>& symmetries, bool verbose);

private:
  static std::vector<int> elements_;
  static bool initialized_;
};

template <class Parameters>
void MomentumExchangeDomain::initialize(const Parameters& parameters, bool verbose) {
  using KDmn = typename Parameters::KClusterDmn::parameter_type;
  using SymmetryDmn = typename domains::cluster_symmetry<KDmn>;
  const auto& symmetry_table = SymmetryDmn::get_symmetry_matrix();
  assert(KDmn::get_size() == symmetry_table[0]);

  // TODO: check for multibands.
  std::vector<std::array<int, 2>> symmetries;
  for (int k1 = 0; k1 < KDmn::get_size(); ++k1)
    for (int sim = 0; sim < symmetry_table[2]; ++sim) {
      const int k2 = symmetry_table(k1, 0, sim).first;
      if (k1 != k2)
        symmetries.push_back(std::array<int, 2>{k1, k2});
    }

  initialize(parameters.compute_all_transfers(), parameters.get_four_point_momentum_transfer_index(),
             KDmn::get_elements(), symmetries, verbose);
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
