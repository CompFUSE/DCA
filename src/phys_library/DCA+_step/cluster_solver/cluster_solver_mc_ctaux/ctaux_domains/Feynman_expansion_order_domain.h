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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_H

#include <vector>

namespace DCA {
namespace QMCI {
// DCA::QMCI::

class Feynman_expansion_order_domain {
  const static int MAX_ORDER_SQUARED = 10000;

public:
  typedef int element_type;

public:
  static int get_size();
  static std::vector<int>& get_elements();

private:
  static std::vector<int>& initialize_elements();
};

int Feynman_expansion_order_domain::get_size() {
  const static int size = MAX_ORDER_SQUARED;
  return size;
}

std::vector<int>& Feynman_expansion_order_domain::get_elements() {
  static std::vector<int>& v = initialize_elements();
  return v;
}

std::vector<int>& Feynman_expansion_order_domain::initialize_elements() {
  static std::vector<int> v(get_size());

  for (int i = 0; i < get_size(); i++)
    v[i] = i;

  return v;
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_DOMAINS_FEYNMAN_EXPANSION_ORDER_DOMAIN_H
