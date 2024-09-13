// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Kernels interface for the class TimeCorrelator.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_KERNELS_INTERFACE_HPP

#include "dca/platform/dca_gpu.h"

#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/device_interpolation_data.hpp"
#include "dca/util/type_help.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

  using dca::util::SignType;
  
template <typename Scalar, typename Real>
void computeG0(linalg::MatrixView<Scalar, linalg::GPU>& g0_mat,
               DeviceInterpolationData<Scalar, SignType<Scalar>> g0, const Real* t_l, const int* b_l,
               const int* r_lf, const Real* t_r, const int* b_r, const int* r_r, const dca::linalg::util::GpuStream& stream);

  

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_KERNELS_INTERFACE_HPP
