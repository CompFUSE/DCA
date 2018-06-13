// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the various minimization methods used for the CPE analytic continuation.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

enum MinimizationMethod {
  GRADIENT_METHOD,
  WEIGHTED_GRADIENT_METHOD,
  CONJUGATE_GRADIENT_METHOD,
  LEVMAR_LIBRARY
};

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP
