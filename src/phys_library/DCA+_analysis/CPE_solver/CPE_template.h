// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Empty class template for a CPE analytic continuation.

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_CPE_TEMPLATE_H
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_CPE_TEMPLATE_H

#include "phys_library/DCA+_analysis/CPE_solver/minimization_method.hpp"

namespace DCA {
template <class parameters_type, class basis_function_t, class k_dmn_t, class w_dmn_t,
          MinimizationMethod minimization_method_t>
class continuous_pole_expansion {};
}

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_CPE_TEMPLATE_H
