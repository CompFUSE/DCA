// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
