// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file defines the types of error computation done by the Monte Carlo solver.

#ifndef DCA_PHYS_ERROR_COMPUTATION_TYPE_HPP
#define DCA_PHYS_ERROR_COMPUTATION_TYPE_HPP

#include <string>

namespace dca {
namespace phys {
// dca::phys::

enum class ErrorComputationType { NONE, STANDARD_DEVIATION, JACK_KNIFE };

ErrorComputationType stringToErrorComputationType(const std::string& str);

std::string toString(ErrorComputationType type);

}  // phys
}  // dca

#endif  // DCA_PHYS_ERROR_COMPUTATION_TYPE_HPP
