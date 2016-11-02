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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAIN_DOMAIN_NAMES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAIN_DOMAIN_NAMES_HPP

#include <string>

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

enum COARSEGRAIN_DOMAIN_NAMES {
  ORIGIN,
  Q_FINE,
  K,
  K_PLUS_Q,
  Q_MINUS_K,
  TETRAHEDRON_K,
  TETRAHEDRON_ORIGIN
};

std::string to_str(COARSEGRAIN_DOMAIN_NAMES NAME);

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAIN_DOMAIN_NAMES_HPP
