// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides distribution strategy tags

#ifndef DCA_DIST_TYPE_HPP
#define DCA_DIST_TYPE_HPP

#include <string>

namespace dca {
enum class DistType { NONE, MPI };

DistType stringToDistType(const std::string& name);
std::string toString(DistType type);
}  // namespace dca

#endif  // DCA_DIST_TYPE_HPP
