// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements electron_spin_domain.hpp.

#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::vector<int> electron_band_domain::flavors_ = {};

}  // namespace domains
}  // namespace phys
}  // namespace dca
