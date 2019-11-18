// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Shrink tools algorithms class.
// It is templated on the device type (CPU|GPU).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_ALGORITHMS_SHRINK_TOOLS_ALGORITHMS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_ALGORITHMS_SHRINK_TOOLS_ALGORITHMS_HPP

#include <cassert>
#include <stdexcept>
#include <vector>

#include "dca/linalg/linalg.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

// Empty class template
template <dca::linalg::DeviceType device_t, typename Real>
class SHRINK_TOOLS_ALGORITHMS {};

// Specialization for CPU
#include "shrink_tools_algorithms_cpu.inc"

#ifdef DCA_HAVE_CUDA
// Specialization for GPU
// Uses SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU> and therefore needs to be included after
// shrink_tools_algorithms_cpu.inc.
#include "shrink_tools_algorithms_gpu.inc"
#endif

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_SHRINK_TOOLS_ALGORITHMS_SHRINK_TOOLS_ALGORITHMS_HPP
