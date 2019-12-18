// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the interpolation of \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.
// It is templated on the device type (CPU|GPU).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP

#include <cassert>
#include <utility>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_template.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_kernels.hpp"
#endif

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

// Empty class template
template <dca::linalg::DeviceType device_t, typename parameters_type, typename Real>
class G0_INTERPOLATION {};

// Specialization for CPU
#include "g0_interpolation_cpu.inc"

#ifdef DCA_HAVE_CUDA
// Specialization for GPU
#include "g0_interpolation_gpu.inc"
#endif

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP
