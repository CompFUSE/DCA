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
// This file defines the names of coarse-graining domains and provides a function to convert them
// into strings.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_NAMES_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_NAMES_HPP

#include <stdexcept>
#include <string>

namespace DCA {

enum COARSEGRAIN_DOMAIN_NAMES {
  ORIGIN,
  Q_FINE,
  K,
  K_PLUS_Q,
  Q_MINUS_K,
  TETRAHEDRON_K,
  TETRAHEDRON_ORIGIN
};

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

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_NAMES_HPP
