// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the different types of vertex frequency domains.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_NAME_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_NAME_HPP

#include <string>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

enum VERTEX_FREQUENCY_NAME {
  COMPACT,
  EXTENDED,
  COMPACT_POSITIVE,
  EXTENDED_POSITIVE,
  COMPACT_SORTED,
  EXTENDED_SORTED,
  EXTENDED_BOSONIC,
  EXTENDED_FERMIONIC,
  CORE_SORTED,
  HIGH_FREQUENCY_SORTED
};

std::string to_str(VERTEX_FREQUENCY_NAME NAME);

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_NAME_HPP
