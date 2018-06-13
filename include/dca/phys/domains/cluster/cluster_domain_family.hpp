// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a cluster domain family consisting of a real space and a momentum space
// domain.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_HPP

#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
class cluster_domain_family {
public:
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_SHAPE SHAPE = S;

  typedef cluster_domain<scalar_type, D, N, REAL_SPACE, S> r_cluster_type;
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;

  template <typename IOWriter>
  static void write(IOWriter& reader);
};

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
template <typename IOWriter>
void cluster_domain_family<scalar_type, D, N, S>::write(IOWriter& writer) {
  writer.open_group(to_str(N));

  {
    writer.open_group(to_str(MOMENTUM_SPACE));

    writer.execute("basis", k_cluster_type::get_basis_vectors());
    writer.execute("super-basis", k_cluster_type::get_super_basis_vectors());
    writer.execute("elements", k_cluster_type::get_elements());

    writer.close_group();
  }

  {
    writer.open_group(to_str(REAL_SPACE));

    writer.execute("basis", r_cluster_type::get_basis_vectors());
    writer.execute("super-basis", r_cluster_type::get_super_basis_vectors());
    writer.execute("elements", r_cluster_type::get_elements());

    writer.close_group();
  }

  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_HPP
