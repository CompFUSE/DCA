// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
