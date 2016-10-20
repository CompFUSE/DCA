// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the symmetry group levels.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_SYMMETRY_GROUP_LEVEL_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_SYMMETRY_GROUP_LEVEL_HPP

#include <string>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

enum symmetry_group_level { UNIT_CELL, SUPER_CELL };
typedef symmetry_group_level symmetry_group_level_type;

std::string to_str(symmetry_group_level name);

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_SYMMETRY_GROUP_LEVEL_HPP
