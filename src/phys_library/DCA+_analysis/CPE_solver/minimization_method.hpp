// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This file defines the various minimization methods used for the CPE analytic continuation.

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP

enum MinimizationMethod {
  GRADIENT_METHOD,
  WEIGHTED_GRADIENT_METHOD,
  CONJUGATE_GRADIENT_METHOD,
  LEVMAR_LIBRARY
};

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_MINIMIZATION_METHOD_HPP
