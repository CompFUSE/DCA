// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_H

#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/cluster/cluster_typedefs.hpp"

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S>
class cluster_domain_family {
public:
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_SHAPE SHAPE = S;

  typedef cluster_domain<scalar_type, D, N, REAL_SPACE, S> r_cluster_type;
  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;

public:
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

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_FAMILY_H
