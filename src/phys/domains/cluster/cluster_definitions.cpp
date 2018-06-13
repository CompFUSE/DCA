// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements cluster_definitions.hpp.

#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::string to_str(int DIMENSION) {
  return std::to_string(DIMENSION);
}

std::string to_str(CLUSTER_NAMES NAME) {
  switch (NAME) {
    case CLUSTER:
      return "CLUSTER";

    case LATTICE_SP:
      return "LATTICE_SP";

    case LATTICE_TP:
      return "LATTICE_TP";

    case TMP_CLUSTER:
      return "TMP_CLUSTER";

    case VASP_LATTICE:
      return "VASP_LATTICE";

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

std::string to_str(CLUSTER_REPRESENTATION NAME) {
  switch (NAME) {
    case REAL_SPACE:
      return "REAL_SPACE";

    case MOMENTUM_SPACE:
      return "MOMENTUM_SPACE";

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

std::string to_str(CLUSTER_SHAPE NAME) {
  switch (NAME) {
    case BRILLOUIN_ZONE:
      return "BRILLOUIN_ZONE";

    case PARALLELLEPIPEDUM:
      return "PARALLELLEPIPEDUM";

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // domains
}  // phys
}  // dca
