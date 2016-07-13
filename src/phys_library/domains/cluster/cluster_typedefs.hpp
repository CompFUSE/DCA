// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file contains enumerations used to specify the different cluster domains. It also provides
// functions to convert these enumerations into strings.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_TYPEDEFS_HPP
#define PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_TYPEDEFS_HPP

#include <string>
#include <stdexcept>

std::string to_str(int DIMENSION) {
  return std::to_string(DIMENSION);
}

enum CLUSTER_NAMES { CLUSTER, LATTICE_SP, LATTICE_TP, TMP_CLUSTER, VASP_LATTICE };

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

enum CLUSTER_REPRESENTATION { REAL_SPACE, MOMENTUM_SPACE };

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

enum CLUSTER_SHAPE { BRILLOUIN_ZONE, PARALLELLEPIPEDUM };

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

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_TYPEDEFS_HPP
