// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Kernel interface for GPU interpolation.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_KERNELS_INTERFACE_HPP
#ifndef DCA_HAVE_GPU
#error "This file requires GPU."
#endif

#include "dca/platform/dca_gpu.h"

#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/device_interpolation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

template <typename Real>
Real interpolateSlow(Real tau, int linindex, const DeviceInterpolationData<Real>& g0);

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_KERNELS_INTERFACE_HPP
