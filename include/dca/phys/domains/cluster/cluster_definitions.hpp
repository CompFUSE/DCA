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
// This file contains enumerations used to specify the different cluster domains. It also provides
// functions to convert these enumerations into strings.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DEFINITIONS_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DEFINITIONS_HPP

#include <string>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::string to_str(int DIMENSION);

enum CLUSTER_NAMES { CLUSTER, LATTICE_SP, LATTICE_TP, TMP_CLUSTER, VASP_LATTICE };
std::string to_str(CLUSTER_NAMES NAME);

enum CLUSTER_REPRESENTATION { REAL_SPACE, MOMENTUM_SPACE };
std::string to_str(CLUSTER_REPRESENTATION NAME);

enum CLUSTER_SHAPE { BRILLOUIN_ZONE, PARALLELLEPIPEDUM };
std::string to_str(CLUSTER_SHAPE NAME);

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DEFINITIONS_HPP
