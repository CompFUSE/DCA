// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements coarsegrain_domain_names.hpp.

#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegrain_domain_names.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

std::string to_str(COARSEGRAIN_DOMAIN_NAMES NAME) {
  switch (NAME) {
    case ORIGIN:
      return "ORIGIN";

    case Q_FINE:
      return "Q-FINE";

    case K:
      return "K";

    case K_PLUS_Q:
      return "K+Q";

    case Q_MINUS_K:
      return "Q-K";

    case TETRAHEDRON_K:
      return "TETRAHEDRON-K";

    case TETRAHEDRON_ORIGIN:
      return "TETRAHEDRON-ORIGIN";

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return "???";
}

}  // clustermapping
}  // phys
}  // dca
