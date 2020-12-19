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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#ifdef DCA_HAVE_CUDA

#include <cuda.h>

#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/device_interpolation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

template <typename Real>
void buildG0Matrix(linalg::MatrixView<Real, linalg::GPU> G0, const int n_init,
                   const bool right_section, DeviceConfiguration config,
                   DeviceInterpolationData<Real> g0_interp, cudaStream_t stream);

// For testing purposes only: computes on the GPU a single value of the interpolated G0 at the
// desired imaginary time "tau" and sector "linindex". Returns the result.
// This method is dominated by latency in a kernel launch and CPU-GPU communication.
template <typename Real>
Real interpolateSlow(Real tau, int linindex, const DeviceInterpolationData<Real>& g0);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
