// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
